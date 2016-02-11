function posterior = lazy_abc(opts)

if nargin <1
opts.max_samples = 10^3;
opts.batches = 100;
opts.num_particles = 2; 
%opts.num_repeats = 10;
opts.continue_prob = 0.1; %probability of continuing an unlikely calculation in lazy abc
opts.max_jump_length = 20; %average jumps beyond 50 microns are seen as unphysical
opts.is_parallel=1;
opts.ss=1;
opts.alpha=0.1;
opts.save_name = 'v2';
end

addpath ../
R = opts.max_samples/opts.batches; %number of samples in each batch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate_synthetic data
InitOpts = opts;
InitOpts.num_particles = 100;
params = [1.16, 0.8, 0.11, 0.42, 0.84, 0.58, 0.01, 0];
data_summary = create_synthetic_data(params,opts);

theta_store = zeros(opts.max_samples,3);
dist_store = zeros(opts.max_samples,1);
a_star_store = zeros(opts.max_samples,1);

theta_batch = zeros(R,3);
dist_batch = zeros(R,1);
a_star_batch = zeros(R,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:opts.batches
fprintf('Starting the %dth batch... \n', k);

z = (k-1)*R;

parfor j=1:R
[early_stopping, a_star, temp_params] = propose_theta_lazy(params,opts)
if ~early_stopping
	dist = continue_lazy(temp_params,data_summary,opts);
elseif rand(1)<a_star
	fprintf('early stopping applies!!! \n');
	dist = continue_lazy(temp_params,data_summary,opts);
else
        fprintf('early stopping applies!!! \n');
	dist = inf;
end

%record theta value
theta_batch(j,:) = temp_params([1,4,6]);
dist_batch(j,:) = dist;
a_star_batch(j,:) = a_star;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta_store(z+(1:R),:) = theta_batch;
dist_store(z+(1:+R),:) = dist_batch;
a_star_store(z+(1:R),:) = a_star_batch;

fprintf('Batch %d complete out of %d. Saved \n \n ',k, opts.batches);
save(sprintf('lazy_abc%s',opts.save_name),'theta_store','dist_store', 'a_star_store');

end

Y = quantile(dist_store,opts.alpha);
post_ind = (dist_store<Y);
posterior = theta_store(post_ind,:);
fprintf('Posterior Saved \n \n ');
save(sprintf('abc_again%s',opts.save_name),'theta_store','dist_store', 'posterior','Y','opts');

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

function [early_stopping, a_star, params] = propose_theta_lazy(params,opts)
%draw params
rr = rand(3,1);
theta = 10.^[2*rr(1)-1, 2*rr(2)-1, -rr(3)]; %prior uniform on log of parameters
params([1,4,6]) = theta;

av_jump = theta(1)/theta(2); %calculate average jump distance
early_stopping = (av_jump>opts.max_jump_length);
a_star = early_stopping*opts.continue_prob + (1-early_stopping);
end

function [dist] = continue_lazy(params,data_summary,opts)

%generate data
q = create_synthetic_data(params,opts);
dist = my_dist(q,data_summary);
end

