function [abc_theta,abc_weights] = ABC_population_monte_carlo

%created 12/5/15
%last edit 12/5/15
%implements a version of ABC via rejection sampling to fit parameters to
%some 'data' which is first simulated from the model
%Model used is velocity jump process with 2 modes
%see Turner 2012 for ABC population monte carlo tutorial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic; 
rng(22);
close all

%Fake parameters
params.nu1 = 0.4;
params.nu2 = 0.1;
params.lambda = 7.7;
params.omega = 1/6;
params.phi = 0.58;
params.x_0=0.5;  
params.theta_0 = 0; 

%Generate fake data
%calculate appropriate summary statistic - we choose MFPT
mfpt_fake_data = mfpt_calculator(params,100,0,65)

%Choose tolerance sequence
num_generations = 3;
delta = [400,200,100,50,25,12];

%At t=1 for first generation
N=50; %N particles at each generation
%create while loop

%set prior
prior_params = [0.4,0.08,7.7,1/6,0.58,0,0.5];
prior_sigma1 = 0.1;
prior_sigma2 = 0.15;
p_indices = [1,5];

n=0;
abc_theta = zeros(N,2);
abc_weights = zeros(N,1);
fprintf('Generation 1 begins\n');
for i=1:N
    mfpt_candidate_data = mfpt_fake_data+delta(1)+1; %initialise greater than tolerance
while abs(mfpt_fake_data - mfpt_candidate_data)>delta(1)    %compare the summary statistic to that from original data
n=n+1;
rr = rand(1,2);

    %simulate parameters from the prior    
%     candidate.nu1 = 0.4 + 0.1*rr(1,j);
%     candidate.nu2 = 0.08 + 0.02*rr(2,j);
%     candidate.lambda = 7.7; % + 2*randn(1);
%     candidate.omega = 1/6; %and this  %+ 0.5*randn(1);
%     candidate.phi = 0.58; %fix this too %+0.2*randn(1);
%     candidate.theta_0 = 0; %these are fixed
%     candidate.x_0 = 0.5; %these are fixed
par_params = prior_params;
par_params(p_indices) = prior_params(p_indices)+[prior_sigma1,prior_sigma2].*(rr-0.5);

    %simulate data using these model parameters
    %Calculate summary statistic (MFPT)
    mfpt_candidate_data = mfpt_calculator(par_params,20,1,4);
    if n>10^5
        fprintf('Thats more than enough parameters to look at for now.\n')
        break
    end
end    %repeat until N acceptances have been made
abc_theta(i,:) = par_params(p_indices);
abc_weights(i) = 1/N;
end
sigma = 2*var(abc_theta);


figure(1);
plot(abc_theta(:,1),abc_theta(:,2),'o');
hold all
grid on
plot(prior_params(p_indices(1)),prior_params(p_indices(2)),'rx','MarkerSize',12);
set(gca, 'fontsize',14);
xlabel('\nu_1');
ylabel('\nu_2');

%First generation for t=1 is done
%Now loop over generations
for tau=2:num_generations
    %store previous theta
    theta_store = abc_theta;
    par_params = repmat(prior_params,N,1);
    par_params(:,p_indices) = theta_store;
    n=0;
    tau
for i=1:N
    mfpt_candidate_data = mfpt_fake_data+delta(tau)+1; %initialise greater than tolerance
while abs(mfpt_fake_data - mfpt_candidate_data)>delta(tau)    %compare the summary statistic to that from original data
n=n+1;

    %simulate parameters from the prior    
%     candidate.nu1 = 0.4 + 0.1*rr(1,j);
%     candidate.nu2 = 0.08 + 0.02*rr(2,j);
%     candidate.lambda = 7.7; % + 2*randn(1);
%     candidate.omega = 1/6; %and this  %+ 0.5*randn(1);
%     candidate.phi = 0.58; %fix this too %+0.2*randn(1);
%     candidate.theta_0 = 0; %these are fixed
%     candidate.x_0 = 0.5; %these are fixed

%sample params from previous iteration
u = rand(1);
my_index = 1;
%draw from discrete distribution with weights abc_weights
while cumsum(abc_weights(1:my_index))/sum(abc_weights)<u
    my_index = my_index+1;
end


%peturb previous parameters
par_params(my_index,p_indices) = par_params(my_index,p_indices) + sigma.*randn(1,2); 

    %simulate data using these model parameters
    %Calculate summary statistic (MFPT)
    mfpt_candidate_data = mfpt_calculator(par_params(my_index,:),20,1,4);
    if n>10^5
        fprintf('Thats more than enough parameters to look at for now.\n')
        break
    end
end    %repeat until N acceptances have been made
abc_theta(i,:) = par_params(my_index,p_indices);
prior = exp(-((abc_theta(i,1)-prior_params(p_indices(1))^2)/(2*prior_sigma1^2)))*exp(-((abc_theta(i,2)-prior_params(p_indices(2)))^2)/(2*prior_sigma2^2));
abc_weights(i) = prior/(sum(abc_weights.*exp(-((abc_theta(:,1)-theta_store(:,1)).^2)/(2*sigma(1)^2))...
    .*exp(-((abc_theta(:,2)-theta_store(:,2)).^2)/(2*sigma(2)^2)))/(sigma(1)*sigma(2))^2);
end
  sigma = 2*var(abc_theta);  
    


abc_theta

figure(1);
plot(abc_theta(:,1),abc_theta(:,2),'o');
hold all
grid on
plot(prior_params(p_indices(1)),prior_params(p_indices(2)),'rx','MarkerSize',12);
set(gca, 'fontsize',14);
xlabel('\nu_1');
ylabel('\nu_2');
end

figure(2);
plot(abc_theta(:,1),abc_theta(:,2),'bo');
hold all
grid on
plot(prior_params(p_indices(1)),prior_params(p_indices(2)),'rx','MarkerSize',12);
set(gca, 'fontsize',14);
xlabel('\nu_1');
ylabel('\nu_2');

toc
end

function mean_fp_time = mfpt_calculator(par_params,num_particles,is_parallel,r_random)
%runs velocityjump2D_ModesInput which is a velocity jump process for a single particle
%runs this many times for many particles and returns simulated estimates of
%mfpt

if is_parallel
    params.nu1 = par_params(1);
    params.nu2 = par_params(2);
    params.lambda = par_params(3);
    params.omega = par_params(4);
    params.phi = par_params(5); 
    params.theta_0 = par_params(6);     
    params.x_0 = par_params(7); 
else
    params = par_params;
end

t_max = 3;

anchor_times = zeros(num_particles,1);
parfor j=1:num_particles
    [~, anchor_times(j), ~,~] = velocityjump2D_ModesInput(t_max, params, 1, 0, 2, 0, r_random*j);
end
mean_fp_time = mean(anchor_times);
%    sd_anchor_time = std(anchor_times);

end
