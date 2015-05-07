function [accepted_params, close_params] = ABC_rejection_sampling

%created 7/5/15
%last edit 7/5/15
%implements a version of ABC via rejection sampling to fit parameters to
%some 'data' which is first simulated from the model
%Model used is velocity jump process with 2 modes
%see Beaumont et al 2002 for ABC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic; 

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
mfpt_fake_data = mfpt_calculator(params,1000,0)

%Choose tolerance
delta = 100;

%create while loop until we find enough good parameter choices

iterations = 100;
n=0;
k=1;
any_accepted = 0;
while any_accepted<k
n=n+1
accepted_params = zeros(iterations,1);
rr = randn(iterations,1);
rs = randn(iterations,1);


parfor j=1:iterations
    %simulate parameters from the prior
    
%     candidate.nu1 = 0.4 + 0.1*rr(1,j);
%     candidate.nu2 = 0.08 + 0.02*rr(2,j);
%     candidate.lambda = 7.7; % + 2*randn(1);
%     candidate.omega = 1/6; %and this  %+ 0.5*randn(1);
%     candidate.phi = 0.58; %fix this too %+0.2*randn(1);
%     candidate.theta_0 = 0; %these are fixed
%     candidate.x_0 = 0.5; %these are fixed

par_params = [0.4+0.1*rr(j),0.08+0.02*rs(j),7.7,1/6,0.58,0,0.5];

    %simulate data using these model parameters
    %Calculate summary statistic (MFPT)
    mfpt_candidate_data = mfpt_calculator(par_params,20,1);
    
    %compare the summary statistic to that from original data
    if abs(mfpt_fake_data - mfpt_candidate_data) < delta
        %if close then accept
        accepted_params(j) = 1;
    end
    %repeat until k acceptances have been made

end
any_accepted = sum(accepted_params);
any_accepted

close_params = [0.4*ones(iterations,1)+0.1*rr,0.08*ones(iterations,1)+0.02*rs, ...
    7.7*ones(iterations,1),1/6*ones(iterations,1), ...
    0.58*ones(iterations,1),0*ones(iterations,1),0.5*ones(iterations,1)].*repmat(accepted_params,1,7);
    
    if n>100
        fprintf('Thats more than enough parameters to look at for now.')
        break
    end
end

close all
figure;
plot(close_params(:,1),close_params(:,2),'bo');
hold on
grid on
plot(params.nu1,params.nu2,'rx');
set(gca, 'fontsize',14);
xlabel('\nu_1');
ylabel('\nu_2');
toc
end

function mean_fp_time = mfpt_calculator(par_params,num_particles,is_parallel)
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

t_max = 2;

anchor_times = zeros(num_particles,1);
for j=1:num_particles
    [~, anchor_times(j), ~,~] = velocityjump2D_ModesInput(t_max, params, 1, 0, 2, 0, 4*j);
end
mean_fp_time = mean(anchor_times);
%    sd_anchor_time = std(anchor_times);

end
