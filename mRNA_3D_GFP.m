function [sum_stat,t_distn,new_coords] = mRNA_3D_GFP(coords, params)
%code wrapped such that it takes list of coords as input and outputs a new
%list of co-ords
%Created 23 Jun 2016
%Last edit 23 Jun 2016
%3D velocity jump model for transport of an RNP along a MT by molecular motors
%based on 2D version
%edited boundary condition so that particles fall off microtubules, rather
%than reflecting
%Assumes: 1 mode of motion, diffusion is neglected
%Now try to use biologically motivated parameters
%Edited to transitions between states are clearer and nucleus is included
%Edited so that parameter input is much nicer
%with_anchoring is logical to include anhoring at site of localization or reflection there
%my_seed is the random seed to be used
%num_modes is for the number of modes of motion ie only AT or with
%diffusion also

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    %fprintf('default parameters used\n')
    params.with_anchoring = 1;
    params.num_modes = 2;
    params.nu1 = 1.16; %speed of RNP complex under active transport [zimyanin et al 2008]
    params.nu2 = 0.80; %ratio between speed for active transport vs diffusion [zimyanin et al 2008]
    params.lambda_1=0;   %1/0.13; %transition rate =7.69 [zimyanin et al 2008]
    params.lambda_2 = 0.11;
    params.omega_1= 0.42;    %1/6*(num_modes>1); %rate of falling off the microtubule [zimyanin et al 2008] since average track length 2.4 - 2.8 microns -> average jump for 6s -> rate 1/6
    params.omega_2 = 0.84;
    params.phi = 0.58; %percentage of microtubules in posterior direction for biased angle distn [parton et al 2011]
    params.Lx = 52; %length of cell in x direction
    params.Ly = 37; %in y direction
    params.Lz = 37; %in z direction
    params.nuc_radius = 10; %radius of nucleus
    params.theta_0 = 0; %initial angle is 0
    params.rc_width = 1;
    params.ellipsoid_boundary = 1;
    params.sample_times = 0:(1.5*10^-4):(4.5*10^-3); %roughly the same as time series for thresholded image
    params.num_hist_bins=21;
    params.zslice = [-1,1];
    params.path_summary_stat=1;
    
    xi = randn(3,780); %78 ish particles in first frame of thresholded image
    lmb = sqrt(sum(xi.^2,1));
    coords = params.nuc_radius*xi./repmat(lmb,3,1); %generate random points at uniform on the unit sphere and scale. Now centred on origin
end
new_coords=zeros(3,size(coords,2));
t_distn = timeseries(zeros(numel(params.sample_times),size(coords,2)),params.sample_times*60^2);
path_stats = zeros(size(coords,2),6);
to_remove=[];
%loop over all particles in the list
for pp=1:size(coords,2)
    is_anchored=0; %has RNP reach destination and anchored yet?
    is_attached=(rand(1)<params.omega_1/(params.omega_1+params.omega_2)); %is molecular motor attached to the microtubule?
    v = params.nu1*is_attached + params.nu2*(1-is_attached); %initialise the speed
    T = params.lambda_2; %store lambda
    params.lambda_2 = (params.num_modes>1)*T; %initially is attached so 0
    nuc_rad = params.nuc_radius^2;
    rc_width = params.rc_width; %3microns wide for ring canals
    time=0; %unit seconds
    
    %select position of next particle in the list
    xpos = coords(:,pp);
    A = diag([1/(params.Lx/2)^2,1/(params.Ly/2)^2,1/(params.Lz/2)^2]);
    %check coords input
    if xpos'*xpos<nuc_rad
        %particle starting in nucleus
        %warning('Please provide particle co-ordinates that do not lie inside the nucleus \n');
    elseif params.ellipsoid_boundary && xpos'*A*xpos>1
        %particle outside ellipsoid domain
        error('Particle outside ellipsoid domain\n');
    elseif ~params.ellipsoid_boundary && (abs(xpos(1))>params.Lx/2 || abs(xpos(2))>params.Ly/2 || abs(xpos(3))>params.Lz/2)
        % particle outside rectangular domain
        error('Particle outside rectangular domain\n');
    end
    
    estimated_max_path_length = round(20000*max(params.sample_times));
    path = zeros(1+estimated_max_path_length,3);
    states = zeros(1+estimated_max_path_length,1);
    states(1) = is_attached;
    path(1,:) = xpos;
    transitions = zeros(estimated_max_path_length,4); %falling off/reattaching events
    jump_times = zeros(1 + estimated_max_path_length,1); %jumps
    n_transition = 0;
    n_jump = 1;
    
    endtime=3600*max(params.sample_times); %one hour.
    %loop through time
    while time<endtime && ~is_anchored
        alpha_rate = params.lambda_1 + params.lambda_2*(1-is_attached) + params.omega_1*is_attached + params.omega_2*(1-is_attached);
        %find next jump
        rr=rand(5,1);
        tau = 1/alpha_rate*log(1/rr(1));
        delx = tau*v;
        phi = acos(2*rr(5)-1);  %cos(phi) must be uniform on [-1,1]
        theta = ((rr(2)<=params.phi)*(rr(2)<=params.phi/2)*rr(3)*pi/2 + ...
            (rr(2)<=params.phi)*(rr(2)>params.phi/2)*(3/2*pi+rr(3)*pi/2) + (rr(2)>params.phi)*(pi/2+pi*rr(3)))*is_attached ...
            + (1-is_attached)*(rr(3))*2*pi;
        
        %update time
        time = time+tau;
        %jumps
        xpos = xpos + delx*[cos(theta)*sin(phi); sin(theta)*sin(phi); cos(phi)]; %delx times unit vec in spherical polars
        
        %check if absorbed at ring canal
        if xpos(1)>=params.Lx/2 && abs(xpos(2))<rc_width && abs(xpos(3))<rc_width %should probably use crossing position, not final position, but effect small if speed smallish
            %anchoring_time = (params.Lx/2-(xpos(1)-delx*cos(theta)*sin(phi)))/(v*cos(theta)*sin(phi)) + time;
            xpos(1) = (params.Lx/2)*(params.with_anchoring>0.5);
            is_anchored = params.with_anchoring; %absorb at right hand boundary or reflect depending on whether we have anchoring
            break; %otherwise may still be reflected
        end
        %check if interacted with boundary
        if params.ellipsoid_boundary==1
            %reflect from an ellipsoid
            if xpos'*A*xpos >1
                %then has interacted with boundary and should fall off MT
                prev_pos = xpos - delx*[cos(theta)*sin(phi); sin(theta)*sin(phi); cos(phi)];
                [xpos,t_intersect] = stop_at_boundary_intersection(prev_pos,xpos,[cos(theta)*sin(phi); sin(theta)*sin(phi); cos(phi)],A);
                time = time - tau + t_intersect/v;
            end
        else
            %assume rectangular/cuboid boundary
            for ii=1:3
                Li = params.Lx*(ii==1) + params.Ly*(ii==2) + params.Lz*(ii==3); %length for x, y and z respectively
                if xpos(ii) >= Li/2   %rectangle
                    direction=[cos(theta)*sin(phi); sin(theta)*sin(phi); cos(phi)];
                    t_extra = (xpos(ii)-Li/2)/direction(ii);
                    time = time - t_extra/v;
                    xpos = xpos - direction*v*t_extra;
                    %                xpos(ii) = Li/2;   %at boundary
                elseif xpos(ii) < -Li/2    %(rectangle)
                    direction=[cos(theta)*sin(phi); sin(theta)*sin(phi); cos(phi)];
                    t_extra = (xpos(ii)+Li/2)/direction(ii);
                    time = time - t_extra/v;
                    xpos = xpos - direction*t_extra;
                    %                xpos(ii) = -Li/2;  %at boundary
                end
            end
        end
        
        %test for interaction with nucleus
        rad = sum(xpos.^2); %euclidean distance from centre of nucleus squared
        if rad<nuc_rad
            %particle has gone inside nucleus so reflect
            prev_pos = xpos - delx*[cos(theta)*sin(phi); sin(theta)*sin(phi); cos(phi)];
            %solve for intersection time
            t_intersect = solve_ellipsoid_intersection(prev_pos,[cos(theta)*sin(phi); sin(theta)*sin(phi); cos(phi)],1/params.nuc_radius^2*eye(3));
            t_intersect = min(t_intersect); %want the first point it intersects
            pos_intersect = prev_pos + t_intersect*[cos(theta)*sin(phi); sin(theta)*sin(phi); cos(phi)];
            xpos = pos_intersect;
            time = time - tau + t_intersect/v;
            rad_check = sum(pos_intersect.^2) - params.nuc_radius^2;
            if abs(rad_check)>10^(-6)
                warning('interacted with nucleus incorrectly'); %no longer an error, as for abc wierd parameter sets can cause bad behaviour, but hopefully these parameters should be discounted
                break; %try to exit time loop to prevent long unphysical simulations
            end
        end
        if rr(4)<params.lambda_2*(1-is_attached)/alpha_rate
            %takes a normal step and changes direction. This has been taken out
            %of the if statement
        elseif rr(4)<(params.lambda_2*(1-is_attached)+params.lambda_1)/alpha_rate
            error('should still have \lambda_1 as 0. Something is up.');
            
        elseif rr(4)<(params.lambda_2*(1-is_attached)+params.lambda_1+params.omega_1*is_attached)/alpha_rate
            if (params.num_modes>1)
                %two modes
                %rates switch as RNP falls off/reattaches onto microtubule
                % falls off
                v = params.nu2;
                is_attached = 0;
                n_transition = n_transition+1;
                transitions(n_transition,:) =  [xpos', time];
            else
                %only one mode
                %take a normal step and change direction. NB outside of if
            end
        elseif rr(4)<(params.lambda_2*(1-is_attached)+params.lambda_1+params.omega_1*is_attached+params.omega_2*(1-is_attached))/alpha_rate
            if (params.num_modes>1)
                %two modes
                %rates switch as RNP falls off/reattaches onto microtubule
                %reattaches
                v = params.nu1;
                is_attached = 1;
                n_transition = n_transition+1;
                transitions(n_transition,:) =  [xpos', time];
            else
                %only one mode
                %take a normal step and change direction. NB outside of if
            end
        else
            %oops
            error('something is wrong')
        end
        
        if time > jump_times(n_jump) + 10^-10 %if we have taken some kind of successful step. Tolerance of 10^-10 to prevent updating when an unsuccessful step is taken on the boundary
            n_jump = n_jump+1;
            path(n_jump,:) = xpos;
            states(n_jump) = is_attached;
            jump_times(n_jump) = time;
        end
    end
    %final outputs for statistics
    
    % if is_anchored
    %     final_position = xpos;
    % else
    %     final_position = xpos;
    %     anchoring_time = inf;
    % end
    %
    % path = path(1:n_jump,:);
    % jump_times = jump_times(1:n_jump);
    if is_anchored
        %add pt at absorbtion
        n_jump = n_jump+1;
        path(n_jump,:) = xpos;
        states(n_jump) = is_attached;
        jump_times(n_jump) = time;
        %and then at end time
        %so don't need to extrapolate for last pt of time series
        n_jump = n_jump+1;
        path(n_jump,:) = xpos;
        states(n_jump) = is_attached;
        jump_times(n_jump) = endtime;
    end
    if params.path_summary_stat
        %calculate path based summary stats
        % be aware of differences between 3D and 2D observations
        %speeds
        framewise_speeds = sqrt(sum(diff(path(1:n_jump,1:3)).^2,2))./diff(jump_times(1:n_jump));
        nu_1 = sum(framewise_speeds.*states(1:n_jump-1))/sum(states(1:n_jump-1));
        nu_2 = sum(framewise_speeds.*(~states(1:n_jump-1)))/sum(~states(1:n_jump-1));
        %mt bias
        framewise_angles = acos(diff(path(1:n_jump,1))./sum(diff(path(1:n_jump,1:3)).^2,2));
        angle_ind = strfind(states(1:n_jump)',[1,0]); %find transitions into AT
        angle = framewise_angles(angle_ind);
        angle = angle(~isnan(angle)); %remove static frames where there is no new angle
        towards_posterior = (angle>pi/2).*(angle<3*pi/2);
        phi = 1 - mean(towards_posterior); %mt bias
        %rate parameters
        k1 = numel(angle_ind);
        k2 = numel(strfind(states(1:n_jump)',[0, 1]));
        k3 = numel(strfind(states(1:n_jump)',[0, 0]));
        k4 = numel(strfind(states(1:n_jump)',[1, 1]));
        %use time spent in each state
        q1 = sum(diff(jump_times(1:n_jump)).*(states(1:n_jump-1)))/sum(states(1:n_jump-1));
        q2 = sum(diff(jump_times(1:n_jump)).*(~states(1:n_jump-1)))/sum(~states(1:n_jump-1));
        z1 = q1*(k4+k1)/k1;
        z2 = q2*(k2+k3)/k2;
        z3 = q2*(k2+k3)/k3;
        if ~isnan([nu_1, nu_2, phi, 1/z1, 1/z2, 1/z3])
            path_stats(pp,:) = [nu_1, nu_2, phi, 1/z1, 1/z2, 1/z3];
        else
            to_remove=[to_remove,pp];
        end
    else
        %use population based distribution for summary stats
        %convert to time series format
        ts = timeseries(path(1:n_jump,:),jump_times(1:n_jump),'Name','mRNA_path');
        %ts1 = resample(ts,params.sample_times*60^2);
        ktemp=resample(ts,params.sample_times*60^2); %resample to discrete times/regular time grid
        zslice_ind = logical((ktemp.Data(:,3)>params.zslice(1))&(ktemp.Data(:,3)<params.zslice(2))); %pick only a single z slice
        t_distn.Data(zslice_ind,pp) = ktemp.Data(zslice_ind,1); %just the xcoord, adjust this if full 3D summary needed
        t_distn.Data(~zslice_ind,pp) = -10*params.Lx; %default value for missing data
    end
    new_coords(:,pp) = xpos;
end %end loop over particles
if n_jump>20000
    warning('Consider increasing estimated path length for efficiency \n');
end
if params.path_summary_stat
    path_stats(to_remove,:)=[];
    sum_stat = mean(path_stats,1);
else
    sum_stat = histc(t_distn.Data',linspace(-params.Lx/2,params.Lx/2,params.num_hist_bins),1);
end
end
