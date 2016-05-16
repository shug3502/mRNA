function flag = CheckaR(a,R)
%a is initial probability vector
%R is transition rate matrix
%%%%%%%%%%%%%%%%%%%%%%%%
if any(a)<0
    fprintf('Oops: some initial state is negative\n');
    flag=1;
elseif abs(sum(a)-1)>10^(-12)
    fprintf('Oops: initial state probs do not sum to 1\n');
    flag=0;
elseif size(R,1)~=size(R,2)
    fprintf('Oops: transition matrix not square\n');
    flag=1;
elseif any(diag(R)>0)
    fprintf('Oops: transition rates must be >=0\n');
    flag=1;
elseif any(sum(R,2)~=0)
    fprintf('Oops: row sums must = 0\n');
    flag=1;
elseif size(R,1)~=length(a)
    fprintf('Oops: rate matrix must match dimensions of inital distn\n');
    flag=1;
else
    flag=0;
end
