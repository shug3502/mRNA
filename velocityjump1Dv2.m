function velocityjump1Dv2
%Created 20 April 2015
%Last edit 20 April 2015
%Simple 1D model for transport of an RNP along a MT by molecular motors
%Assumes: 2 modes of motion, motor can reverse

%isn't this more of a position jump process?

close all
clear all
rng(5); % set random seed

dx=0.01; %compartment width
v = 0.03; %speed of RNP complex under active transport
sigmaV = 20; %ratio between speed for active transport vs diffusion

TR=2; %transition rate right
sigmaR=2; %ratio between transition rates for diffusion/active transport rates
TL=0.5; %transition rate left
sigmaL=.5; %ratio between transition rates for diffusion/active transport rates
switching_rate=0.1; %rate of falling off the microtubule

x_0=dx/2;  %Initially in first compartment, ie. at NPC
is_anchored=0; %has RNP reach destination and anchored yet?
is_attached=1; %is molecular motor attached to the microtubule?

time=0; %unit seconds
xpos = x_0; %initialise x position
path = x_0; %add first point to path
transitions = []; %falling off/reattaching events
jump_times = 0; %jumps 
endtime=100;
%loop through time
while time<endtime && ~is_anchored
    alpha = TR+TL+switching_rate;
    %find next jump
    rr=rand(2,1);
    tau = 1/alpha*log(1/rr(1));
    delx = tau*v;
    
    if rr(2)<TR/alpha
        %jumps right
        xpos = xpos + delx;
        if xpos >=1
            is_anchored = 1; %absorb at right hand boundary
        end 
    elseif rr(2)<(TR+TL)/alpha
        %jumps left
        xpos = xpos - delx;
        if xpos < 0
            xpos = -xpos; %reflect at left hand boundary
        end
    elseif rr(2)<(TR+TL+switching_rate)/alpha
        %rates switch as RNP falls off/reattaches onto microtubule
        if is_attached
            % falls off
            TR = TR/sigmaR;
            TL = TL/sigmaL;
            v = v/sigmaV;
            is_attached = 0;
            transitions = [transitions; xpos, time];
        else
            %reattaches
            TR = TR*sigmaR;
            TL = TL*sigmaL;
            v = v*sigmaV;
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

figure;
plot(jump_times, path, 'r', 'linewidth',3)
xlabel('t')
ylabel('x position along microtubule')
grid on
hold on
plot(transitions(:,2),transitions(:,1),'bo') 
%axis([0 50 0 1])
end



