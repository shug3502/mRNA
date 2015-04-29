function [is_anchored, anchoring_time, final_position_x] = velocityjump1D_fixed(input_time,phi_input,plot_option, random_seed)
%Created 26 April 2015
%Last edit 26 April 2015
%Simple 1D model for transport of an RNP along a MT by molecular motors
%Assumes: 1 modes of motion, motor can reverse


%close all
tic;
rng(random_seed); % set random seed

v_AT = 0.4; %speed of RNP complex under active transport
v_diff = 0.08; %speed of RNP diffusing
v=v_AT; %initially on a MT under active transport

gamma = 1/0.13; % rate of falling off and reattaching
phi = phi_input; % percentage bias in microtubule orientation
switching_rate = 1/6;
domain_size = 80; %80 microns

x_0=0.5;  %Initially at NPC
is_anchored=0; %has RNP reach destination and anchored yet?
is_attached=1; %is molecular motor attached to the microtubule?

time=0; %unit seconds
xpos = x_0; %initialise x position
path = x_0; %add first point to path
transitions = []; %falling off/reattaching events
jump_times = 0; %jumps 
endtime=input_time*60^2; %input time in hrs, endtime in seconds
%loop through time
while time<endtime && ~is_anchored
    alpha = gamma+switching_rate;
    %find next jump
    rr=rand(3,1);
    tau = 1/alpha*log(1/rr(1));
    theta = (rr(2)<=phi/100)*rand(1)*pi + (rr(2)>phi/100)*(pi+pi*rand(1));
    delx = tau*v*sin(theta);
    
    if rr(3)<gamma/alpha
        %jumps right
        xpos = xpos + delx;
        if xpos >=domain_size
            is_anchored = 1; %absorb at right hand boundary
        elseif xpos < 0
            xpos = -xpos; %reflect at left hand boundary
        end
    elseif rr(3)<(gamma+switching_rate)/alpha
        %rates switch as RNP falls off/reattaches onto microtubule
        if is_attached
            % falls off
            v=v_diff;
            is_attached = 0;
            transitions = [transitions; xpos, time];
        else
            %reattaches
            v=v_AT;
            is_attached = 1;
            transitions = [transitions; xpos, time];
        end
    else
        %oops
        error('something is wrong')
    end
    time = time+tau;
    
    path = [path, xpos];
    jump_times = [jump_times, time];
end
anchoring_time = time;
final_position_x = xpos;

if plot_option
    if is_anchored
        plot_col = 'r';
    else
        plot_col = 'g';
    end
    figure(1);
    plot(jump_times, path, plot_col, 'linewidth',3)
    set(gca,'fontsize',16)
    xlabel('t')
    ylabel('x position along microtubule')
    grid on
    hold on
    
    if ~isempty(transitions)
        plot(transitions(:,2),transitions(:,1),'bo')
    end
end
toc
end



