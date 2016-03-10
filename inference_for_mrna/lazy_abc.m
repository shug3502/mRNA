function posterior = lazy_abc(opts)

if nargin <1
opts.max_samples = 10^5;
opts.batches = 100;
opts.num_particles = 5; 
%opts.num_repeats = 10;
opts.continue_prob = 0.1; %probability of continuing an unlikely calculation in lazy abc
opts.max_jump_length = 20; %average jumps beyond 50 microns are seen as unphysical
opts.is_parallel=1;
opts.ss=0;
opts.alpha=0.01;
opts.save_name = 'combined_v2';
end

addpath ../
R = opts.max_samples/opts.batches; %number of samples in each batch
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate_synthetic data
InitOpts = opts;
InitOpts.num_particles = 50;
params = [1.16, 0.8, 0.11, 0.42, 0.84, 0.58, 0.01, 0];
data_summary = create_synthetic_data(params,InitOpts);

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
save(sprintf('lazy_abc%s',opts.save_name),'theta_store','dist_store', 'posterior','Y','opts');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

