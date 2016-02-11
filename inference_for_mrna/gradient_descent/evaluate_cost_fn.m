function cost = evaluate_cost_fn(theta,data_summary,opts)

%created 1/2/2016 JH
%last edit 2/2/2016

%assume sum of squares for now

%addpath ../../
%defaults for calculating summary stat:
% for par_params expect a vector of params of length 7
is_parallel=1;
%option_a_summary_statistic=1; %use distribution at various times as summary stat

p_indices = [1,4,6];
params = [1.16, 0.8, 0.11, 0.42, 0.84, 0.58, 0.01, 0];
params(p_indices) = exp(theta);

%sz1=53; %from size of q_estimate
%sz2=21; %from size of q_estimate
%fprintf('Calculating summary statistics ...\n');

c=zeros(opts.num_repeats,1);
for i=1:opts.num_repeats
    %generate summary stat
    q = summary_statistic_calculator_3D(params,opts.num_particles,is_parallel,opts.ss);
    M = reshape(q,1,[]);
c(i) = norm(M - data_summary);
end
cost = mean(c);
end
