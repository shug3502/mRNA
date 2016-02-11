%script to produce a file with values of theta used to generate certain
%data (summary stats of some data) from the 3D mRNA model
%Edited to use the 3D model, simplified with one mode
%Last edit 12/1/16

%setup
addpath ../
num_samples = 10^4;

%defaults for calculating summary stat:
% for par_params expect a vector of params of length 7
num_particles=100; %only one repeat of data needed in calculating summary stat for a given theta?
is_parallel=1;
option_a_summary_statistic=0; %use distribution at various times as summary stat

%set uniform prior
prior_params = [1, 0, 0, 1, 0, 0.5, 0.01];
p_indices = [1,4,6];
prior_sigma = [2, 0, 0, 2, 0, 1, 0.02];
% prior_params = [1.16, 0.8, 0.11, 0.42, 0.84, 0.58, 0.01, 0];
par_params = prior_params;
% par_params = prior_params(p_indices);

sz1=53; %from size of q_estimate
sz2=21; %from size of q_estimate
M = zeros(num_samples,sz1*sz2+length(p_indices));
fprintf('Calculating summary statistics ...\n');
for j=1:num_samples
	if mod(j,100)==0
		fprintf('%f percent done \n',j/num_samples*100);
	end
    %sample theta from the prior
    rr = rand(1,length(p_indices));
    par_params(p_indices) = prior_params(p_indices)+prior_sigma(p_indices).*(rr-0.5);
    
    %generate summary stat
    q = summary_statistic_calculator_3D(par_params,num_particles,is_parallel,option_a_summary_statistic);
    M(j,1:length(p_indices)) = par_params(p_indices);
    M(j,(1+length(p_indices):end)) = reshape(q,1,[]);
end
fprintf('All done!\n');
csvwrite('train3D100particles.csv',M);

