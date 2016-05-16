%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x=samplefromp(p,n)
%Inputs - p is the probability vector of length k
% - n is the number of random
% integers from 1,2, ...,k returned
%Output - a row vector of length n with entries
% from the set {1, 2, ..., k} with
% probabilities specified by p.
%source http://www.math.wustl.edu/~feres/Math450Lect01.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(p,1)~=1
    error('expected a row vector');
end
k=size(p,2);
u=rand(1,n);
x=zeros(1,n);
for i=1:k
x=x+i*(sum(p(1:i))-p(i) <=u & u<sum(p(1:i)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%