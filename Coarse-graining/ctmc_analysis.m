function [distn,event_times,y] = ctmc_analysis(initial_distn,end_time)
%created 10/5/16
%last edit 11/5/16
%JH
%ctmc with 16 states

%initial_distn should be a row vector, t is time
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(initial_distn,1)~=1
    error('initial_distn should be a row vector \n');
end

%define parameters
params.dt = 0.1;
params.production_rate=0.1;
params.flow_rate=1;
params.save_name = 'coarse-grained-NCs-ctmc-model';

num_time_steps = end_time/params.dt;
%%%%%%%%%%%%%%%
tic;
can_exit_oocyte=0;
Q = get_nc_transitions(can_exit_oocyte); %get rate matrix
Q = break_connections(Q,1,1); %break connections of given NC
if CheckaR(initial_distn,Q); %check this
    error('Error in initial distn or Q matrix\n')
end
producers = eye(16); %producers(1,1)=0;
Q = params.flow_rate*Q; %+params.production_rate*producers;

%%%%%%%%%%%%%%%%%%
%update distn via matrix exponentiation
P = expm(Q*params.dt); %probability of transitions in dt
distn=zeros(num_time_steps+1,size(initial_distn,2));
distn(1,:) = initial_distn;
for tau=1:num_time_steps
    distn(tau+1,:) = distn(tau,:)*P+params.production_rate*params.dt*diag(producers)';
end
toc
if params.production_rate==0
%simulate (one path of) the process directly
[event_times, y]=ctmc_feres(end_time,initial_distn,Q);
toc;
end
distn=distn(1:10:end,:);
save(params.save_name);

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

function Q = break_connections(Q,index1,index2)
%pick a cell to block a given RC both IN and OUT
Q(index1,index2)=0;
Q(index2,index1)=0;
r = 1:16; r(index1)=[];
Q(index1,index1)=-sum(Q(index1,r));
r = 1:16; r(index2)=[];
Q(index2,index2)=-sum(Q(index2,r));








