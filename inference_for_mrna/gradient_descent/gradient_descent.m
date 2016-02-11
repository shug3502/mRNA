function output = gradient_descent(opts)

if nargin <1
opts.theta0 = [0,0,-log(2)];
opts.max_steps = 10^3;
opts.gamma = 0.9;
opts.alpha = 10^-5; %step size
opts.delta = 10^-5;
opts.num_particles = 2; 
opts.num_repeats = 10;
opts.ss=1;
opts.save_name = 'v1';
end

addpath ../../

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initialise
theta = opts.theta0;
step = 0;
is_converged=0; %boolean for if optimization has converged
cost = zeros(1,opts.max_steps);
j=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generate_synthetic data
params = [1.16, 0.8, 0.11, 0.42, 0.84, 0.58, 0.01, 0];
data_summary = create_synthetic_data(params,opts);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while ~is_converged && j<opts.max_steps

%evaluate
cost(j) = evaluate_cost_fn(theta,data_summary,opts);

%estimate gradient
gradient_est = estimate_gradient(theta,data_summary,opts,cost(j));

%update with momentum
step = opts.gamma*step + opts.alpha*gradient_est;
theta = theta - step; %descent step
j=j+1; %update number of steps taken
fprintf('Step %d, theta (%f, %f, %f) \n',j,exp(theta(1)),exp(theta(2)),exp(theta(3)));

if mod(j,100)==0
fprintf('Steps taken: %d \n',j);
fprintf('Saved \n \n ');
save(sprintf('gradient_descent%s',opts.save_name),'cost');
opts.alpha = opts.alpha/2;
end

if norm(step)<opts.delta %then has converged?
is_converged = 1;
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output = cost;
end

function data_summary = create_synthetic_data(params,opts)

%defaults for calculating summary stat:
% for par_params expect a vector of params of length 7
is_parallel=1;
option_a_summary_statistic=1; %use distribution at various times as summary stat
    %generate summary stat
    q = summary_statistic_calculator_3D(params,opts.num_particles,is_parallel,opts.ss);
    data_summary = reshape(q,1,[]);
end
