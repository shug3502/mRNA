function [is_anchored, anchoring_time, final_position_x, final_position_y, pathx, jump_times] = velocityjump2D_ModesInput(input_time, params, with_anchoring, absorb_in_region, num_modes, with_plot, my_seed)
%Created 20 April 2015
%Last edit 7 May 2015
%2D velocity jump model for transport of an RNP along a MT by molecular motors
%Assumes: 1 mode of motion, diffusion is neglected
%Now try to use biologically motivated parameters

%Edited so that parameter input is much nicer

%input_time is length of stage 9 ie 6 to 10 in units of hours
%with_anchoring is logical to include anhoring at site of localization or reflection there
%with_plot is logical to include plotting or not
%my_seed is the random seed to be used
%num_modes is for the number of modes of motion ie only AT or with
%diffusion also
if nargin ~= 7
    fprintf('default parameters used\n')
    input_time = 1;
    with_anchoring = 1;
    absorb_in_region = 0;
    num_modes = 1;
    with_plot = 1;
    my_seed = 32;
    %Now with data on nurse cells from Alex Davidson
    params.nu1 = 1.16; %speed of RNP complex under active transport [zimyanin et al 2008]
    params.nu2 = 0.80; %ratio between speed for active transport vs diffusion [zimyanin et al 2008]
    params.lambda=0.42;   %1/0.13; %transition rate =7.69 [zimyanin et al 2008]
    params.omega=1/6*(num_modes>1); %rate of falling off the microtubule [zimyanin et al 2008] since average track length 2.4 - 2.8 microns -> average jump for 6s -> rate 1/6
    params.phi = 0.58; %percentage of microtubules in posterior direction for biased angle distn [parton et al 2011]
    params.x_0=0.5;  %Initially in first compartment, ie. at NPC
    params.theta_0 = 0; %initial angle is 0
end
%rng(my_seed); % set random seed

%dx=1; %compartment width for initial condition

domain_size = 30; %30 microns for a nurse cell (R Parton) %80 microns [zimyanin et al 2008]
a = domain_size/2;
b=a/2;

y_0 = 0;  %domain_size/4; %initially half way up y axis
is_anchored=0; %has RNP reach destination and anchored yet?
is_attached=1; %is molecular motor attached to the microtubule?
v = params.nu1*is_attached + params.nu2*(1-is_attached); %initialise the speed
T = params.lambda; %store lambda
params.lambda = (1-is_attached*(num_modes>1))*T; %initially is attached so 0

time=0; %unit seconds
xpos = params.x_0; %initialise x position
ypos = y_0; %initialise y

%pathx = [params.x_0; zeros(10^6,1)]; %add first point to path
%pathy = [y_0; zeros(10^6,1)];
estimated_max_path_length = round(50000*input_time);
pathx = zeros(1+estimated_max_path_length,1);
pathx(1) = params.x_0;
pathy = zeros(1+estimated_max_path_length,1);
pathy(1) = y_0;
transitions = zeros(estimated_max_path_length,3); %falling off/reattaching events
jump_times = zeros(1 + estimated_max_path_length,1); %jumps
n_transition = 0;
n_jump = 1;

endtime=3600*input_time; %one hour. Note that stage 9 of development lasts 6-10 hrs - or is this 6-10 since start? [zimyanin et al 2008]
%loop through time
while time<endtime && ~is_anchored
    alpha = params.lambda + params.omega*(num_modes>1);
    %find next jump
    rr=rand(3,1);
    tau = 1/alpha*log(1/rr(1));
    delx = tau*v;
    theta = ((rr(2)<=params.phi)*rand(1)*pi + (rr(2)>params.phi)*(pi+pi*rand(1)))*is_attached ...
        + (1-is_attached)*(rand(1))*2*pi;
    
    if time == 0
        theta = params.theta_0; %set an initial angle. Get rid of this to have T(theta) initially
    end
    %jumps
    xpos = xpos + delx*sin(theta);
    ypos = ypos + delx*cos(theta);
    if xpos >=domain_size-absorb_in_region*2
        anchoring_time = (domain_size-(xpos-delx*sin(theta)))/(v*sin(theta))+ time;
        xpos = (2*domain_size - xpos)*(with_anchoring<0.5)+(domain_size-absorb_in_region*2)*(with_anchoring>0.5);
        is_anchored = 1*with_anchoring; %absorb at right hand boundary or reflect depending on whether we have anchoring
    elseif xpos < 0
        xpos = -xpos; %reflect at left hand boundary
    end
    if ypos >= b*sqrt(1-((xpos-domain_size/2)/a)^2)  %domain_size/2 %rectangle
        ypos = 2*b*sqrt(1-((xpos-domain_size/2)/a)^2) - ypos;  %domain_size - ypos; %reflect at boundary
    elseif ypos < -b*sqrt(1-((xpos-domain_size/2)/a)^2)   % -domain_size/2; %(rectangle)
        ypos = -2*b*sqrt(1-((xpos-domain_size/2)/a)^2)-ypos; %-domain_size - ypos; %reflect at boundary
    end
    
    if rr(3)<params.lambda/alpha
        %takes a normal step and changes direction. This has been taken out
        %of the if statement
    elseif rr(3)<(params.lambda + params.omega*(num_modes>1))/alpha
        %rates switch as RNP falls off/reattaches onto microtubule
        if is_attached
            % falls off
            v = params.nu2;
            is_attached = 0;
            params.lambda = T; %switching within phase when in diffusion
            %if with_plot
            n_transition = n_transition+1;
            transitions(n_transition,:) =  [xpos, ypos, time];
            %end
        else
            %reattaches
            v = params.nu1;
            is_attached = 1;
            params.lambda = T*(num_modes<2); % no switching within phase when in AT
            %if with_plot
            n_transition = n_transition+1;
            transitions(n_transition,:) =  [xpos, ypos, time];
            %end
        end
    else
        %oops
        error('something is wrong')
    end
    time = time+tau;
    %If needing to plot thun uncomment these lines!!
    %    if with_plot
    n_jump = n_jump+1;
    pathx(n_jump) = xpos;
    pathy(n_jump) = ypos;
    jump_times(n_jump) = time;
    %    end
end
%final outputs for statistics
if is_anchored
    final_position_x = domain_size;
else
    final_position_x = xpos;
    anchoring_time = inf;
end

pathx = pathx(1:n_jump);
pathy = pathy(1:n_jump);
jump_times = jump_times(1:n_jump);
final_position_y = ypos;
if with_plot
    %close all
    %some plots
    
    if is_anchored == 1
        plot_col = 'r--'; %plot in red if reached target within simulation time
    else
        plot_col = 'g--'; %green otherwise
    end
    figure(1);
    plot3(jump_times, pathx, pathy, plot_col,'linewidth',3)
    set(gca,'fontsize',14)
    xlabel('t');
    ylabel('x position along microtubule');
    zlabel('y position');
    grid on
    hold all
    if ~isempty(transitions)
        plot3(transitions(:,3),transitions(:,1),transitions(:,2),'bo')
    end
    %plot3(transitions(:,3),transitions(:,1),b*sqrt(1-((transitions(:,1)-domain_size/2)/a).^2),'g')
    %plot3(transitions(:,3),transitions(:,1),-b*sqrt(1-((transitions(:,1)-domain_size/2)/a).^2),'g')
    plot3(0:3600/(domain_size*2):3600,0:0.5:domain_size,b*sqrt(1-(((0:0.5:domain_size)-domain_size/2)/a).^2),'k-')
    plot3(0:3600/(domain_size*2):3600,0:0.5:domain_size,-b*sqrt(1-(((0:0.5:domain_size)-domain_size/2)/a).^2),'k-')
    
    figure(3)
    plot(jump_times, pathx, plot_col, 'linewidth',3)
    set(gca,'fontsize',16)
    xlabel('t');
    ylabel('x position along microtubule');
    grid on
    hold all
    if ~ isempty(transitions)
        plot(transitions(:,3),transitions(:,1),'bo')
    end
    mu_initial = 0.5; %initial mean position
    nu_1 = 0.4; %speed in active transport mode
    F1 = (4*params.phi-2)/pi;
    t = jump_times;
    mu1 = min((nu_1*F1*t + mu_initial),domain_size); %note time scaled to seconds
    if num_modes>1
        mu1 = min((0.5*nu_1*F1*t + mu_initial),domain_size); %adjust for different analytical results for different numbers of modes
    end
    plot(t, mu1, 'b--', 'linewidth',3);
    set(gca, 'fontsize',14)
end

%toc
end



