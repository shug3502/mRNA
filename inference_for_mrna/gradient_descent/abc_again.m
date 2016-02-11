function posterior = abc_again(opts)

if nargin <1
opts.max_samples = 10^3;
opts.num_particles = 2; 
%opts.num_repeats = 10;
opts.is_parallel=1;
opts.ss=1;
opts.alpha=0.1;
opts.save_name = 'v1';
end

addpath ../../

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialise
j=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate_synthetic data
InitOpts = opts;
InitOpts.num_particles = 100;
params = [1.16, 0.8, 0.11, 0.42, 0.84, 0.58, 0.01, 0];
data_summary = create_synthetic_data(params,InitOpts);

theta_store = zeros(opts.max_samples,3);
dist_store = zeros(opts.max_samples,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
parfor j=1:opts.max_samples

[dist, theta] = propose_theta(params,data_summary,opts);

%record theta value
theta_store(j,:) = theta;
dist_store(j,:) = dist;


if mod(j,100)==0
fprintf('Samples generated: %d \n',j);
fprintf('Saved \n \n ');
%save(sprintf('abc_again%s',opts.save_name),'theta_store','dist_store');
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y = quantile(dist_store,opts.alpha);
post_ind = (dist_store<Y);
posterior = theta_store(post_ind,:);

end

function data_summary = create_synthetic_data(params,opts)

%defaults for calculating summary stat:
% for par_params expect a vector of params of length 7
%is_parallel=1;
%option_a_summary_statistic=1; %use distribution at various times as summary stat
    %generate summary stat
    q = summary_statistic_calculator_3D(params,opts.num_particles,opts.is_parallel,opts.ss);
    data_summary = reshape(q,1,[]);
end

function D = my_dist(data1,data2)
scale=1;
D = norm((data1-data2)./scale,2);
end

function [dist, theta] = propose_theta(params,data_summary,opts)
%draw params
rr = rand(3,1);
theta = [2*rr(1)-1, 2*rr(2)-1, -rr(3)];
params([1,4,6]) = theta;

%generate data
q = create_synthetic_data(params,opts);
dist = my_dist(q,data_summary);
end
