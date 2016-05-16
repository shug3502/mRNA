%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t y]=ctmc_feres(end_time,pi,Q)
%Obtain a sample path with n events for a
%continuous-times Markov chain with initial
%distribution pi and generator matrix Q.
%The output consists of two row vectors:
%the event times t and the vector of states y.
%Vectors t and y may be shorter than n if
%an absorbing state is found before event n.
%Uses samplefromp(pi,n).
%source http://www.math.wustl.edu/~feres/Math450Lect05.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t=0;
y=samplefromp(pi,1); %initial state
k=0;
while t(k+1)<end_time
    k=k+1; %keep track of number of steps
    i=y(k);
    q=-Q(i,i);
    if q==0
        break
    else
        s=-log(rand)/(-Q(i,i)); %exponential holding time
        t=[t t(k)+s];
        p=Q(i,:);
        p(i)=0;
        p=p/sum(p);
        y=[y samplefromp(p,1)];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%