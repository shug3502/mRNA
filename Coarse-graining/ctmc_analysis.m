function [distn,event_times,y] = ctmc_analysis(initial_distn,params)
%created 10/5/16
%last edit 11/5/16
%JH
%ctmc with 16 states

%initial_distn should be a row vector, t is time
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(initial_distn,1)~=1
    error('initial_distn should be a row vector \n');
end
if isempty(params)
%define (default) parameters
params.end_time = 10;
params.dt = 0.1;
params.production_rate=0;
params.flow_rate=1;
params.thin_results = 10;
params.save_name = 'coarse-grained-NCs-ctmc-model';
params.store_output = 1;
params.connections_to_break = [1, 2]; %ie don't break anything
end

num_time_steps = params.end_time/params.dt;
%%%%%%%%%%%%%%%
can_exit_oocyte=0;
Q = get_nc_transitions(can_exit_oocyte); %get rate matrix
Q = break_connections(Q,params.connections_to_break(1),params.connections_to_break(2)); %break connections of given NC
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
%%%%%%%%%%%%%%%%%%
%can simulate directly
if params.production_rate==0
    %simulate (one path of) the process directly
    [event_times, y]=ctmc_feres(params.end_time,initial_distn,Q);
    toc;
end
distn=distn(1:params.thin_results:end,:);
if params.store_output
    save(params.save_name);
end


