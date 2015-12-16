%script to produce a file with values of theta used to generate certain
%data (summary stats of some data) from the 2D mRNA model
%Last edit 16/12/15

%setup
addpath ../
num_samples = 10^5;

%defaults for calculating summary stat:
% for par_params expect a vector of params of length 7
num_particles=10; %only one repeat of data needed in calculating summary stat for a given theta?
is_parallel=1;
option_a_summary_statistic=0; %use distribution at various times as summary stat

%set uniform prior
prior_params = [1, 1, 1, 1, 1, 0.5, 0.01, 0];
p_indices = 1:7;
prior_sigma = [1, 1, 1, 1, 1, 1, 0.02, 0.02];
% prior_params = [1.16, 0.8, 0.42, 0.42, 0.84, 0.58, 0.01, 0];
par_params = prior_params(p_indices);

sz1=53; %from size of q_estimate
sz2=21; %from size of q_estimate
M = zeros(num_samples,sz1*sz2+length(p_indices));

for j=1:num_samples
    %sample theta from the prior
    rr = rand(1,length(p_indices));
    par_params(p_indices) = prior_params(p_indices)+prior_sigma(p_indices).*(rr-0.5);
    
    %generate summary stat
    [q_estimate] = summary_statistic_calculator(par_params,num_particles,is_parallel,option_a_summary_statistic);
    M(j,1:length(p_indices)) = par_params(p_indices);
    M(j,(1+length(p_indices):end)) = reshape(q_estimate,1,[]);
end

csvwrite('train.csv',M);
