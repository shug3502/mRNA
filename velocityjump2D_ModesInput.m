function [is_anchored, anchoring_time, final_position_x, final_position_y] = velocityjump2D_ModesInput(input_time, phi_input, x_initial_input, theta_0_input, with_anchoring, absorb_in_region, num_modes, with_plot, my_seed)
%Created 20 April 2015
%Last edit 20 April 2015
%2D velocity jump model for transport of an RNP along a MT by molecular motors
%Assumes: 1 mode of motion, diffusion is neglected
%Now try to use biologically motivated parameters

%input_time is length of stage 9 ie 6 to 10 in units of hours
%with_anchoring is logical to include anhoring at site of localization or reflection there
%with_plot is logical to include plotting or not
%my_seed is the random seed to be used
%num_modes is for the number of modes of motion ie only AT or with
%diffusion also

rng(my_seed); % set random seed

%tic;
dx=1; %compartment width for initial condition
domain_size = 30; %30 microns for a nurse cell (R Parton) %80 microns [zimyanin et al 2008]
a = domain_size/2;
b=a/2;

v = 0.4; %speed of RNP complex under active transport [zimyanin et al 2008]
sigmaV = 5; %ratio between speed for active transport vs diffusion [zimyanin et al 2008]

T=1/0.13; %transition rate =7.69 [zimyanin et al 2008]
switching_rate=1/6*(num_modes>1); %rate of falling off the microtubule [zimyanin et al 2008] since average track length 2.4 - 2.8 microns -> average jump for 6s -> rate 1/6 
phi1 = phi_input; %58; %percentage of microtubules in posterior direction for biased angle distn [parton et al 2011]

x_0=x_initial_input;  %dx/2;  %Initially in first compartment, ie. at NPC
y_0 = 0;  %domain_size/4; %initially half way up y axis
is_anchored=0; %has RNP reach destination and anchored yet?
is_attached=1; %is molecular motor attached to the microtubule?

time=0; %unit seconds
xpos = x_0; %initialise x position
ypos = y_0; %initialise y
pathx = x_0; %add first point to path
pathy = y_0;
transitions = []; %falling off/reattaching events
jump_times = 0; %jumps 
endtime=3600*input_time; %one hour. Note that stage 9 of development lasts 6-10 hrs - or is this 6-10 since start? [zimyanin et al 2008]
%loop through time
while time<endtime && ~is_anchored
    alpha = T+switching_rate;
    %find next jump
    rr=rand(3,1);
    tau = 1/alpha*log(1/rr(1));
    delx = tau*v;
    theta = ((rr(2)<=phi1/100)*rand(1)*pi + (rr(2)>phi1/100)*(pi+pi*rand(1)))*is_attached ... 
        + (1-is_attached)*(rand(1))*2*pi;
    
    if time == 0
        theta = theta_0_input; %set an initial angle. Get rid of this to have T(theta) initially
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
    
    if rr(3)<T/alpha
        %takes a normal step and changes direction. This has been taken out
        %of the if statement
    elseif rr(3)<(T+switching_rate)/alpha
        %rates switch as RNP falls off/reattaches onto microtubule
        if is_attached
            % falls off
            v = v/sigmaV;
            is_attached = 0;
            transitions = [transitions; xpos, ypos, time];
        else
            %reattaches
            v = v*sigmaV;
            is_attached = 1;
            transitions = [transitions; xpos, ypos, time];
        end
    else
        %oops
        error('something is wrong')
    end
    time = time+tau;
    
    pathx = [pathx, xpos];
    pathy = [pathy, ypos];
    jump_times = [jump_times, time];
end
%final outputs for statistics
if is_anchored
    final_position_x = domain_size;
else
    final_position_x = xpos;
    anchoring_time = inf;
end

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
end

%toc
end



