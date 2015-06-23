function [posterior] = update_frame_to_frame(state, position, angle, num_repeats, params)

dt=1/3*60^(-2); %a third of a second between frames
    %set parameters
    params.nu1 = 1.16; %speed of RNP complex under active transport [zimyanin et al 2008]
    params.nu2 = 0.80; %ratio between speed for active transport vs diffusion [zimyanin et al 2008]
    params.lambda_1=0;   %1/0.13; %transition rate =7.69 [zimyanin et al 2008]
    params.lambda_2 = 0.11;
    params.omega_1= 0.42;    %1/6*(num_modes>1); %rate of falling off the microtubule [zimyanin et al 2008] since average track length 2.4 - 2.8 microns -> average jump for 6s -> rate 1/6
    params.omega_2 = 0.84;
    params.phi = 0.58; %percentage of microtubules in posterior direction for biased angle distn [parton et al 2011]
    params.x_0=position(1);  %Initially in first compartment, ie. at NPC
    params.y_0=position(2);
    params.Lx = 52; %length of cell in x direction
    params.Ly = 37; %in y direction
    params.nuc_radius = 10; %radius of nucleus
    params.theta_0 = angle; %initial angle is 0

for i=1:num_repeats
[~, ~, final_position_x, final_position_y, pathx, pathy, jump_times, is_attached, theta] = aux_velocityjump2D(dt, params, 1, state, 2)

end

end

function [is_anchored, anchoring_time, final_position_x, final_position_y, pathx, pathy, jump_times, is_attached, theta] = aux_velocityjump2D(input_time, params, with_anchoring, state, num_modes)
%Created 20 April 2015
%Last edit 22 May 2015
%Based on velocityjump2D_with_nucleus
%2D velocity jump model for transport of an RNP along a MT by molecular motors
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
if nargin ~= 5
    fprintf('default parameters used\n')
    input_time = 1;
    with_anchoring = 1;
    num_modes = 1;
    %Now with data on nurse cells from Alex Davidson
    params.nu1 = 1.16; %speed of RNP complex under active transport [zimyanin et al 2008]
    params.nu2 = 0.80; %ratio between speed for active transport vs diffusion [zimyanin et al 2008]
    params.lambda_1=0;   %1/0.13; %transition rate =7.69 [zimyanin et al 2008]
    params.lambda_2 = 0.11;
    params.omega_1= 0.42;    %1/6*(num_modes>1); %rate of falling off the microtubule [zimyanin et al 2008] since average track length 2.4 - 2.8 microns -> average jump for 6s -> rate 1/6
    params.omega_2 = 0.84;
    params.phi = 0.58; %percentage of microtubules in posterior direction for biased angle distn [parton et al 2011]
    params.x_0=0.5;  %Initially in first compartment, ie. at NPC
    params.y_0=0.5;
    params.Lx = 52; %length of cell in x direction
    params.Ly = 37; %in y direction
    params.nuc_radius = 10; %radius of nucleus
    params.theta_0 = 0; %initial angle is 0
end
%rng(3); % set random seed

is_anchored=0; %has RNP reach destination and anchored yet?
is_attached=state; %Put this equal to previous state value  %is molecular motor attached to the microtubule?
v = params.nu1*is_attached + params.nu2*(1-is_attached); %initialise the speed
T = params.lambda_2; %store lambda
params.lambda_2 = (1-is_attached*(num_modes>1))*T; %initially is attached so 0
nuc_rad = params.nuc_radius^2;
rc_width = 1; %2microns wide for ring canals

time=0; %unit seconds
xpos =  params.x_0;  %params.Lx/2 + params.nuc_radius*(2*rand(1)-1); %initialise x position
ypos = params.y_0;   %sqrt(params.nuc_radius^2-(xpos-params.Lx/2).^2); %initialise y

%pathx = [params.x_0; zeros(10^6,1)]; %add first point to path
%pathy = [y_0; zeros(10^6,1)];
estimated_max_path_length = round(20000*input_time);
pathx = zeros(1+estimated_max_path_length,1);
pathx(1) = xpos;
pathy = zeros(1+estimated_max_path_length,1);
pathy(1) = ypos;
transitions = zeros(estimated_max_path_length,3); %falling off/reattaching events
jump_times = zeros(1 + estimated_max_path_length,1); %jumps
n_transition = 0;
n_jump = 1;

endtime=3600*input_time; %one hour.
%loop through time
while time<endtime && ~is_anchored
    alpha = params.lambda_1 + params.lambda_2*(1-is_attached) + params.omega_1 + params.omega_2;
    %find next jump
    rr=rand(4,1);
    tau = 1/alpha*log(1/rr(1));
    delx = tau*v;
    theta = ((rr(2)<=params.phi)*rr(3)*pi + (rr(2)>params.phi)*(pi+pi*rr(3)))*is_attached ...
        + (1-is_attached)*(rr(3))*2*pi;
    
        if time == 0
            theta = params.theta_0; %set an initial angle. Get rid of this to have T(theta) initially
        end
    
    %jumps
    xpos = xpos + delx*sin(theta);
    ypos = ypos + delx*cos(theta);
    
    if xpos >=params.Lx
        if abs(ypos)<rc_width %should probably use crossing position, not final position, but effect small if speed smallish
        anchoring_time = (params.Lx-(xpos-delx*sin(theta)))/(v*sin(theta))+ time;
        xpos = (2*params.Lx - xpos)*(with_anchoring<0.5)+(params.Lx)*(with_anchoring>0.5);
        is_anchored = with_anchoring; %absorb at right hand boundary or reflect depending on whether we have anchoring
        else
        xpos = 2*params.Lx - xpos;    
        if xpos<0
            xpos = min(-xpos,params.Lx);
        end
        end
    elseif xpos < 0
        xpos = -xpos; %reflect at left hand boundary
    end
    if ypos >= params.Ly/2  %rectangle
        ypos = params.Ly - ypos;   %reflect at boundary
    elseif ypos < -params.Ly/2    %(rectangle)
        ypos = -params.Ly-ypos;  %reflect at boundary
    end
    
    %test for interaction with nucleus
    rad = (xpos-params.Lx/2)^2 + ypos^2; %euclidean distance from centre of nucleus squared
    if rad<nuc_rad
        %particle has gone inside nucleus so reflect
        G = 1/tan(theta); %gradient
        x_intersection = [((params.Lx - 2*G*(ypos-G*xpos)) -sqrt((2*G*(ypos-G*xpos)-params.Lx)^2 - 4*(1+G^2)*((params.Lx^2)/4+(ypos-G*xpos)^2-params.nuc_radius^2)))/(2*(1+G^2)),...
            ((params.Lx - 2*G*(ypos-G*xpos))+sqrt((2*G*(ypos-G*xpos)-params.Lx)^2 - 4*(1+G^2)*((params.Lx^2)/4+(ypos-G*xpos)^2-params.nuc_radius^2)))/(2*(1+G^2))];
        [~,index] = min((x_intersection - xpos).^2);
        x_intersection = x_intersection(index);
        y_intersection =  sign(ypos)*sqrt(params.nuc_radius^2-(x_intersection-params.Lx/2).^2);   %ypos + G*(x_intersection-xpos)
        %check
        rad_check = (x_intersection - params.Lx/2)^2 + y_intersection^2 - params.nuc_radius^2;
        if rad_check>10^(-6)
            x1 = [xpos - delx*sin(theta),ypos - delx*cos(theta)]
            x2 = [xpos,ypos]
            error('interacted with nucleus incorrectly');
        end

        reflect_dist = sqrt((ypos-y_intersection)^2+(xpos-x_intersection)^2);
        xpos  = x_intersection + sin(2*pi-theta)*reflect_dist;
        ypos = y_intersection - cos(2*pi-theta)*reflect_dist;
    end
    
    if rr(4)<params.lambda_2*(1-is_attached)/alpha
        %takes a normal step and changes direction. This has been taken out
        %of the if statement
    elseif rr(4)<(params.lambda_2*(1-is_attached)+params.lambda_1)/alpha
        error('should still have \lambda_1 as 0. Something is up.');
        
    elseif rr(4)<(params.lambda_2*(1-is_attached)+params.lambda_1+params.omega_1)/alpha
        if (num_modes>1)
            %two modes
            %rates switch as RNP falls off/reattaches onto microtubule
            % falls off
            v = params.nu2;
            is_attached = 0;
            n_transition = n_transition+1;
            transitions(n_transition,:) =  [xpos, ypos, time];
        else
            %only one mode
            %take a normal step and change direction. NB outside of if
        end
    elseif rr(4)<(params.lambda_2*(1-is_attached)+params.lambda_1+params.omega_1+params.omega_2)/alpha
        if (num_modes>1)
            %two modes
            %rates switch as RNP falls off/reattaches onto microtubule
            %reattaches
            v = params.nu1;
            is_attached = 1;
            n_transition = n_transition+1;
            transitions(n_transition,:) =  [xpos, ypos, time];
        else
            %only one mode
            %take a normal step and change direction. NB outside of if
        end
    else
        %oops
        error('something is wrong')
    end
    time = time+tau;
    n_jump = n_jump+1;
    pathx(n_jump) = xpos;
    pathy(n_jump) = ypos;
    jump_times(n_jump) = time;
end
%final outputs for statistics
if is_anchored
    final_position_x = params.Lx;
else
    final_position_x = xpos;
    anchoring_time = inf;
end

pathx = pathx(1:n_jump);
pathy = pathy(1:n_jump);
jump_times = jump_times(1:n_jump);
final_position_y = ypos;


%toc
end



