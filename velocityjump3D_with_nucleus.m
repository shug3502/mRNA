function [is_anchored, anchoring_time, final_position, path, jump_times] = velocityjump3D_with_nucleus(params)
%Created 4 Jan 2016
%Last edit 4 Jan 2016
%3D velocity jump model for transport of an RNP along a MT by molecular motors
%based on 2D version
%Assumes: 1 mode of motion, diffusion is neglected
%Now try to use biologically motivated parameters

%Edited to transitions between states are clearer and nucleus is included

%Edited so that parameter input is much nicer

%input_time is length of stage 9 ie 6 to 10 in units of hours
%with_anchoring is logical to include anhoring at site of localization or reflection there
%with_plot is logical to include plotting or not
%my_seed is the random seed to be used
%num_modes is for the number of modes of motion ie only AT or with
%diffusion also

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin ~= 1
    %fprintf('default parameters used\n')
    params.input_time = 10;
    params.with_anchoring = 1;
    params.num_modes = 2;
    params.with_plot = 0;
    %Now with data on nurse cells from Alex Davidson
    params.nu1 = 1.16; %speed of RNP complex under active transport [zimyanin et al 2008]
    params.nu2 = 0.80; %ratio between speed for active transport vs diffusion [zimyanin et al 2008]
    params.lambda_1=0;   %1/0.13; %transition rate =7.69 [zimyanin et al 2008]
    params.lambda_2 = 0.11;
    params.omega_1= 0.42;    %1/6*(num_modes>1); %rate of falling off the microtubule [zimyanin et al 2008] since average track length 2.4 - 2.8 microns -> average jump for 6s -> rate 1/6
    params.omega_2 = 0.84;
    params.phi = 0.58; %percentage of microtubules in posterior direction for biased angle distn [parton et al 2011]
    params.x_0=0.5;  %Initially in first compartment, ie. at NPC
    params.Lx = 52; %length of cell in x direction
    params.Ly = 37; %in y direction
    params.Lz = 37; %in z direction
    params.nuc_radius = 10; %radius of nucleus
    params.theta_0 = 0; %initial angle is 0
    params.rc_width = 1;
    params.ellipsoid_boundary = 1;
end
if params.with_plot
    rng(493); % set random seed
end

is_anchored=0; %has RNP reach destination and anchored yet?
is_attached=1; %is molecular motor attached to the microtubule?
v = params.nu1*is_attached + params.nu2*(1-is_attached); %initialise the speed
T = params.lambda_2; %store lambda
params.lambda_2 = (1-is_attached*(params.num_modes>1))*T; %initially is attached so 0
nuc_rad = params.nuc_radius^2;
rc_width = params.rc_width; %3microns wide for ring canals

time=0; %unit seconds
xi = randn(3,1);
lmb = sqrt(sum(xi.^2));
xpos = params.nuc_radius*xi/lmb; %generate random points at uniform on the unit sphere and scale. Now centred on origin.
%xpos = [-params.nuc_radius;0;0];

estimated_max_path_length = round(200000*params.input_time);
path = zeros(1+estimated_max_path_length,3);
path(1,:) = xpos;
transitions = zeros(estimated_max_path_length,4); %falling off/reattaching events
jump_times = zeros(1 + estimated_max_path_length,1); %jumps
n_transition = 0;
n_jump = 1;

endtime=3600*params.input_time; %one hour.
%loop through time
while time<endtime && ~is_anchored
    alpha = params.lambda_1 + params.lambda_2*(1-is_attached) + params.omega_1 + params.omega_2;
    %find next jump
    rr=rand(5,1);
    tau = 1/alpha*log(1/rr(1));
    delx = tau*v;
    phi = acos(2*rr(5)-1);  %cos(phi) must be uniform on [-1,1]
    theta = ((rr(2)<=params.phi)*(rr(2)<=params.phi/2)*rr(3)*pi/2 + ...
        (rr(2)<=params.phi)*(rr(2)>params.phi/2)*(3/2*pi+rr(3)*pi/2) + (rr(2)>params.phi)*(pi/2+pi*rr(3)))*is_attached ...
        + (1-is_attached)*(rr(3))*2*pi;
    
    %     if time == 0
    %         theta = params.theta_0; %set an initial angle. Get rid of this to have T(theta) initially
    %     end
    
    %update time
    time = time+tau;
    %jumps
    xpos = xpos + delx*[cos(theta)*sin(phi); sin(theta)*sin(phi); cos(phi)]; %delx times unit vec in spherical polars
    
    if xpos(1) >=params.Lx/2
        if abs(xpos(2))<rc_width && abs(xpos(3))<rc_width %should probably use crossing position, not final position, but effect small if speed smallish
            anchoring_time = (params.Lx/2-(xpos(1)-delx*cos(theta)*sin(phi)))/(v*cos(theta)*sin(phi)) + time;
            xpos(1) = (params.Lx - xpos(1))*(params.with_anchoring<0.5)+(params.Lx/2)*(params.with_anchoring>0.5);
            is_anchored = params.with_anchoring; %absorb at right hand boundary or reflect depending on whether we have anchoring
        else
            xpos(1) = params.Lx - xpos(1);
            if xpos(1)<-params.Lx/2
                xpos(1) = -params.Lx-xpos(1);
            end
        end
    end
    if params.ellipsoid_boundary==1
        %reflect from an ellipsoid
        A = diag([1/(params.Lx/2)^2,1/(params.Ly/2)^2,1/(params.Lz/2)^2]);
        if xpos'*A*xpos >1
            %then has interacted with boundary and needs reflecting
            prev_pos = xpos - delx*[cos(theta)*sin(phi); sin(theta)*sin(phi); cos(phi)];
            t_intersect = solve_ellipsoid_intersection(prev_pos,[cos(theta)*sin(phi); sin(theta)*sin(phi); cos(phi)],A);
            if isnan(t_intersect)
                error('t is nan');
            end
            t_intersect = 1/v*max(t_intersect); %want the first point it intersects
            pos_intersect = prev_pos + t_intersect*v*[cos(theta)*sin(phi); sin(theta)*sin(phi); cos(phi)];
            unit_normal = -diag(A)'*pos_intersect; %inward pointing normal?
            unit_normal = unit_normal./norm(unit_normal); %normalise
            
            Lstar = 2*pos_intersect - xpos;
            xpos = 2*(unit_normal'*Lstar)*unit_normal - Lstar; %new reflected position
            if xpos'*A*xpos -1 > 10^-6 %some tolerance
                xpos = pos_intersect; %take position as intersection pt. NB Not treated properly for second reflection
            end
        end
    else
        %assume rectangular/cuboid boundary
        for ii=1:3
            Li = params.Lx*(ii==1) + params.Ly*(ii==2) + params.Lz*(ii==3); %length for x, y and z respectively
            if xpos(ii) >= Li/2 && ii~=1   %rectangle
                xpos(ii) = Li - xpos(ii);   %reflect at boundary
            elseif xpos(ii) < -Li/2    %(rectangle)
                xpos(ii) = -Li-xpos(ii);  %reflect at boundary
            end
        end
    end
    %test for interaction with nucleus
    rad = sum(xpos.^2); %euclidean distance from centre of nucleus squared
    if rad<nuc_rad
        %particle has gone inside nucleus so reflect
        
        prev_pos = xpos - delx*[cos(theta)*sin(phi); sin(theta)*sin(phi); cos(phi)];
        %solve for intersection time
        %bb = ((prev_pos(1)*cos(theta)+prev_pos(2)*sin(theta))*sin(phi)+prev_pos(3)*cos(phi));
        %t_intersect = [-bb+sqrt(bb^2-sum(prev_pos.^2)+params.nuc_radius^2); -bb-sqrt(bb^2-sum(prev_pos.^2)+params.nuc_radius^2)]  %scaled
        t_intersect = solve_ellipsoid_intersection(prev_pos,[cos(theta)*sin(phi); sin(theta)*sin(phi); cos(phi)],1/params.nuc_radius^2*eye(3));
        t_intersect = 1/v*min(t_intersect); %want the first point it intersects
        pos_intersect = prev_pos + t_intersect*v*[cos(theta)*sin(phi); sin(theta)*sin(phi); cos(phi)];
        %unit_normal = -pos_intersect/params.nuc_radius;
        unit_normal = params.nuc_radius*ones(1,3)*pos_intersect;
        unit_normal = unit_normal./norm(unit_normal);
        
        Lstar = 2*pos_intersect - xpos;
        xpos = 2*(unit_normal'*Lstar)*unit_normal - Lstar; %new reflected position
        
        rad_check = sum(pos_intersect.^2) - params.nuc_radius^2;
        if abs(rad_check)>10^(-6)
            error('interacted with nucleus incorrectly');
        end
        
        %        hold on
        %        grid on
        %        plot(linspace(-params.nuc_radius,params.nuc_radius,50),sqrt(params.nuc_radius^2-(linspace(params.Lx/2-params.nuc_radius,params.Lx/2+params.nuc_radius,50)-params.Lx/2).^2),'k');
        %        plot(linspace(-params.nuc_radius,+params.nuc_radius,50),-sqrt(params.nuc_radius^2-(linspace(params.Lx/2-params.nuc_radius,params.Lx/2+params.nuc_radius,50)-params.Lx/2).^2),'k');
        %        axis([0,params.Lx,-params.Ly/2,params.Ly/2]);
        %        set(gca,'fontsize',14);
        %        xlabel('x')
        %        ylabel('y')
        %        plot(xpos,ypos,'ro','MarkerSize',12); %red in x_2
        %        plot(x_intersection,y_intersection,'bo','MarkerSize',12); %blue is intersection x*
        %        hold on
        %        plot(xpos,ypos,'go','MarkerSize',12); % green is x_2'
        %        error('enough now');
    end
    
    if rr(4)<params.lambda_2*(1-is_attached)/alpha
        %takes a normal step and changes direction. This has been taken out
        %of the if statement
    elseif rr(4)<(params.lambda_2*(1-is_attached)+params.lambda_1)/alpha
        error('should still have \lambda_1 as 0. Something is up.');
        
    elseif rr(4)<(params.lambda_2*(1-is_attached)+params.lambda_1+params.omega_1)/alpha
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
    elseif rr(4)<(params.lambda_2*(1-is_attached)+params.lambda_1+params.omega_1+params.omega_2)/alpha
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
    
    n_jump = n_jump+1;
    path(n_jump,:) = xpos;
    jump_times(n_jump) = time;
end
%final outputs for statistics
if is_anchored
    final_position = xpos;
else
    final_position = xpos;
    anchoring_time = inf;
end

path = path(1:n_jump,:);
jump_times = jump_times(1:n_jump);
if params.with_plot
    close all
    %some plots
    
    transitions = transitions(1:n_transition,:);
    
    if is_anchored == 1
        plot_col = 'r'; %plot in red if reached target within simulation time
    else
        plot_col = 'g'; %green otherwise
    end
    %     figure(1);
    %     plot3(jump_times, path(:,1), path(:,2), plot_col,'linewidth',3)
    %     set(gca,'fontsize',14)
    %     xlabel('t');
    %     ylabel('x position');
    %     zlabel('y position');
    %     grid on
    %     hold all
    %     if ~isempty(transitions)
    %         plot3(transitions(:,3),transitions(:,1),transitions(:,2),'bo')
    %     end
    %     plot3(linspace(0,0,50),linspace(params.Lx/2-params.nuc_radius,params.Lx/2+params.nuc_radius,50),sqrt(params.nuc_radius^2-(linspace(params.Lx/2-params.nuc_radius,params.Lx/2+params.nuc_radius,50)-params.Lx/2).^2),'k-')
    %     plot3(linspace(0,0,50),linspace(params.Lx/2-params.nuc_radius,params.Lx/2+params.nuc_radius,50),-sqrt(params.nuc_radius^2-(linspace(params.Lx/2-params.nuc_radius,params.Lx/2+params.nuc_radius,50)-params.Lx/2).^2),'k-')
    %     %plot3(transitions(:,3),transitions(:,1),b*sqrt(1-((transitions(:,1)-domain_size/2)/a).^2),'g')
    %     %plot3(transitions(:,3),transitions(:,1),-b*sqrt(1-((transitions(:,1)-domain_size/2)/a).^2),'g')
    %     %     plot3(0:3600/(params.Lx*2):3600,0:0.5:L,b*sqrt(1-(((0:0.5:L)-L/2)/a).^2),'k-')
    %     %     plot3(0:3600/(L*2):3600,0:0.5:L,-b*sqrt(1-(((0:0.5:L)-L/2)/a).^2),'k-')
    
    figure(3)
    subplot(2,1,1)
    plot(jump_times, path(:,1), plot_col, 'linewidth',2)
    set(gca,'fontsize',16)
    xlabel('t');
    ylabel('x position');
    grid on
    hold all
    if ~ isempty(transitions)
        plot(transitions(:,end),transitions(:,1),'ko','markersize',3)
        %axis equal
        axis([0 transitions(end,end) -params.Lx/2 params.Lx/2]);
    end
    set(gca, 'fontsize',18)
    
    %figure(5)
    subplot(2,1,2)
    plot(path(:,1),path(:,2),plot_col,'linewidth',2)
    set(gca,'fontsize',18)
    xlabel('x position');
    ylabel('y position');
    %grid on
    hold all
    if ~ isempty(transitions)
        plot(transitions(:,1),transitions(:,2),'ko','markersize',3)
    end
    hold on
    nuc_vec = linspace(-params.nuc_radius,params.nuc_radius,2000);
    fill(nuc_vec, sqrt(params.nuc_radius^2-(nuc_vec).^2),'k', nuc_vec, -sqrt(params.nuc_radius^2-(nuc_vec).^2),'k');
    
    hold on
    plot((params.Lx/2)*ones(50,1),linspace(-rc_width,rc_width,50),'c','linewidth',5)
    rectangle('Position',[-params.Lx/2,-params.Ly/2,params.Lx,params.Ly],...
        'LineWidth',2,...
        'LineStyle','--')
    axis equal
    %axis([-params.Lx/2, params.Lx/2, -20, 20]);
    axis([-60, 60, -20, 20]);
    %axis equal
    
    print('Figures_for_writeup/Typical_path3D','-depsc');
    %geometry only
    figure(4);
    set(gca,'fontsize',18)
    xlabel('x position');
    ylabel('y position');
    %grid on
    hold on
    nuc_vec = linspace(-params.nuc_radius,params.nuc_radius,2000);
    fill(nuc_vec, sqrt(params.nuc_radius^2-(nuc_vec).^2),'k', nuc_vec, -sqrt(params.nuc_radius^2-(nuc_vec).^2),'k');
    hold on
    plot((params.Lx/2)*ones(50,1),linspace(-rc_width,rc_width,50),'c','linewidth',5)
    rectangle('Position',[-params.Lx/2,-params.Ly/2,params.Lx,params.Ly],...
        'LineWidth',2,...
        'LineStyle','--')
    axis equal
    axis([-60, 60, -20, 20]);
    %axis equal
    print('Figures_for_writeup/Domain_geometry3D','-depsc');
    
    
    %code to make a movie
    %     nframe=min(200,length(path));
    %     mov(1:nframe)=struct('cdata',[],'colormap',[]);
    %     figure;
    %     hold on
    %     grid on
    %     plot(linspace(-params.nuc_radius,params.nuc_radius,50),sqrt(params.nuc_radius^2-(linspace(-params.nuc_radius,params.nuc_radius,50)).^2),'k');
    %     plot(linspace(-params.nuc_radius,params.nuc_radius,50),-sqrt(params.nuc_radius^2-(linspace(-params.nuc_radius,params.nuc_radius,50)).^2),'k');
    %     axis([0,params.Lx,-params.Ly/2,params.Ly/2]);
    %     set(gca,'fontsize',14);
    %     xlabel('x')
    %     ylabel('y')
    %     for i=1:nframe
    %         frame_num = i;
    %         plot(pathx(frame_num),pathy(frame_num),'ro','MarkerSize',12);
    %         mov(frame_num)=getframe(gcf);
    %     end
    %     movie2avi(mov, 'VJ_With_Nucleus.avi', 'compression', 'None');
end

%toc
end



