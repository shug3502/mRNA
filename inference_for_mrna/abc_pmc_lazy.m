function posterior = abc_pmc_lazy(opts)
%population monte carlo with adaptive distance function, see prangle 2015
%last edit 18/2/2016

if nargin <1
    opts.N = 400;
    opts.alpha = 0.5;
    opts.max_generations = 10;
    opts.num_particles = 5;
    opts.continue_prob = 0.1; %probability of continuing an unlikely calculation in lazy abc
    opts.max_jump_length = 20; %average jumps beyond 50 microns are seen as unphysical
    opts.is_parallel=1;
    opts.min_ess = opts.N/4;
    opts.save_name = 'v4';
end

rng(111)
addpath ../
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate_synthetic data
InitOpts = opts;
InitOpts.num_particles = 100;
params = [1.16, 0.8, 0.11, 0.42, 0.84, 0.58, 0.01, 0];
data_summary = create_synthetic_data(params,InitOpts);



theta_store = zeros(opts.N,3,opts.max_generations);
dist_store = zeros(opts.N,opts.max_generations);
a_star_store = zeros(opts.N,opts.max_generations);
ss_store = zeros(opts.N,numel(data_summary),opts.max_generations);
weights_store = zeros(opts.N,opts.max_generations);
h_store = zeros(1,opts.max_generations+1);

R = opts.N/opts.alpha; %upper bound estimated number of samples in each batch
theta_batch = zeros(R,3);
dist_batch = zeros(R,1);
a_star_batch = zeros(R,1);
ss_batch = zeros(R,numel(data_summary));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t = 1; %set generation as 1 initially
h_store(1) = 50; %inf %for first generation, accept all
%M = opts.N;     M = ceil(opts.N/opts.alpha); %number of particles to accept at each generation
ss_scale=1; weights=1; sigma=1;

while t<=opts.max_generations  %or some other condition
    num_accepted = 0;
    accepted = zeros(R,1);
    j = 0;
    
    fprintf('Starting generation %d ... \n',t);
    while num_accepted < opts.N
        j = j+1;
        if mod(j,100)==0
            fprintf('Done %d iterations, number accepted: %d \n',j,num_accepted);
        end
        [early_stopping, a_star, temp_params] = propose_theta_lazy(params,opts,t, theta_store, weights,sigma);
        if ~early_stopping
            q = continue_lazy(temp_params,opts);
            dist = calculate_distance(q,data_summary,ss_scale);
        elseif rand(1)<a_star
            fprintf('early stopping applies!!! \n');
            q = continue_lazy(temp_params,opts);
            dist = calculate_distance(q,data_summary,ss_scale);
        else
            fprintf('early stopping applies!!! \n');
            dist = inf;
            if ~exist('q','var')  %if not yet defined
                q = data_summary;
            end
        end
        %record theta value
        theta_batch(j,:) = temp_params([1,4,6]);
        dist_batch(j,:) = dist;
        a_star_batch(j,:) = a_star;
        ss_batch(j,:) = q;
        accepted(j) = (dist<h_store(t));
        num_accepted = num_accepted+accepted(j);
    end
    
    h_store(t+1) = h_store(t)*opts.alpha    %quantile(dist_batch(1:j),opts.alpha)
    ind  = logical(accepted);      %ind = (dist_batch<h); %find which samples to keep
    theta_store(:,:,t) = theta_batch(ind,:);
    dist_store(:,t) = dist_batch(ind);
    a_star_store(:,t) = a_star_batch(ind);
    ss_store(:,:,t) = ss_batch(ind,:);
    
    %now calculate the weights
    if t==1
        weights = ones(opts.N,1)./a_star_batch(ind);    
        weights_store(:,t) = weights;
    else
        ii = 0;
        weights = ones(opts.N,1);
        pp = prior_probability(theta_batch(ind,:));
        for i=1:j
            if ind(i)
                ii=ii+1;
                kernel = prod(normpdf(repmat(log10(theta_batch(i,:)),opts.N,1),log10(theta_store(:,:,t-1)),repmat(sigma,opts.N,1)),2);
                qq = sum(weights_store(:,t-1).*kernel)./sum(weights_store(:,t-1));
                weights(ii) = pp(ii)./(qq*a_star_batch(ind(i)));
            end
        end
        weights_store(:,t) = weights;
        norm_weights = weights./sum(weights);
        ESS = 1/sum(norm_weights.^2)
        if ESS < opts.min_ess
            warning('Degenerate: effective sample size too small. Resampling ... \n');
            weights_store(:,t) = ones(opts.N,1);
        end
    end

    sigma = 2*var(log10(theta_store(:,:,t)),weights_store(:,t));  %use log of parameters
    ss_scale= std(ss_store(:,:,t)); %scale for summary statistics - this may be 0 for many statistics
    ss_scale(abs(ss_scale)<10^-8) = 1; % if the scale is zero (ie no data), then use a scale of 1
    
    fprintf('Batch %d complete out of %d. Saved \n \n ',t, opts.max_generations);
    save(sprintf('abc_pmc_lazy%s',opts.save_name));
    
    
    t = t+1;  %to the next generation
    
end

t = opts.max_generations;
Y = quantile(dist_store(:,t),opts.alpha);
post_ind = (dist_store(:,t)<Y);
posterior = theta_store(post_ind,:,t);
fprintf('Run complete \n \n ');
save(sprintf('abc_pmc_lazy%s',opts.save_name));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data_summary = create_synthetic_data(params,opts)

%defaults for calculating summary stat:
% for par_params expect a vector of params of length 71:opts.N
%is_parallel=1;
%option_a_summary_statistic=1; %use distribution at various times as summary stat
%generate summary stat
q = summary_statistic_calculator_combined_3D(params,opts.num_particles,opts.is_parallel);
%   q = summary_statistic_calculator_3D(params,opts.num_particles,opts.is_parallel,opts.ss);
data_summary = reshape(q,1,[]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function D = my_dist(data1,data2,scale)
scale = 1;
D = norm((data1-data2)./scale,2);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [early_stopping, a_star, params] = propose_theta_lazy(params,opts,t, theta_store, weights,sigma)
%draw params
if t==1
    %use prior
    rr = rand(3,1);
    theta = 10.^[2*rr(1)-1, 2*rr(2)-1, -rr(3)]; %prior uniform on log of parameters
else
    in_prior = 0;
    while ~in_prior  % do again until you pick something in the support of the prior
        %perturb previous params
        sample_ind = discretesample(weights,1);
        theta = 10.^(log10(theta_store(sample_ind,:,t-1)) + randn(1,3).*sigma); %sigma should be empirical varianbce of previous generation
        %check if still in prior
        in_prior = min((log10(theta)>[-1,-1,-1]).*(log10(theta)<[1,1,0]));
        %that is nu and lambda between 0.1 and 10, phi between 0.1 and 1
    end
    
end
params([1,4,6]) = theta;
av_jump = theta(1)/theta(2); %calculate average jump distance
early_stopping = (av_jump>opts.max_jump_length);
a_star = early_stopping*opts.continue_prob + (1-early_stopping);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function q = continue_lazy(params,opts)

%generate data
q = create_synthetic_data(params,opts);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dist = calculate_distance(q,data_summary,scale)
dist = my_dist(q,data_summary,scale);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function p = prior_probability(theta)
p = prod(1./(2*theta),2);
end
