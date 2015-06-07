function [abc_theta,abc_weights,entropy,quality,contained_in_pred_interval] = Sample_From_Prior(N,option_a,my_seed)

%created 4/6/15
%last edit 28/5/15
%AS simple as ABC
%Based on Simplest_ABC_with_moments
%added in knn selection of theta for parameters rather than less than a
%specific delta value
%implements a version of ABC via population mc to fit parameters to
%some 'data' which is first simulated from the model
%Model used is velocity jump process with 1 mode
%see Turner 2012 for ABC population monte carlo tutorial
%Should be extendable to adaptive popMc-ABC
%See Lenormand 2013 for APMC
%Should be edited to now use new code for vel jump with nucleus

%Simple rejection sampling ABC on complicated model to match other code for
%APMC etc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
rng(my_seed);
close all

%Fake parameters
params.nu1 = 1.0; %speed of RNP complex under active transport [zimyanin et al 2008]
params.nu2 = 0.80; %ratio between speed for active transport vs diffusion [zimyanin et al 2008]
params.lambda_1=0;   %1/0.13; %transition rate =7.69 [zimyanin et al 2008]
params.lambda_2 = 0.11;
params.omega_1= 0.42;    %1/6*(num_modes>1); %rate of falling off the microtubule [zimyanin et al 2008] since average track length 2.4 - 2.8 microns -> average jump for 6s -> rate 1/6
params.omega_2 = 0.84;
params.phi = 0.58; %0.58 %percentage of microtubules in posterior direction for biased angle distn [parton et al 2011]
params.x_0=0.5;  %Initially in first compartment, ie. at NPC
params.Lx = 52; %length of cell in x direction
params.Ly = 37; %in y direction
params.nuc_radius = 10; %radius of nucleus
params.theta_0 = 0; %initial angle is 0

real_params = [params.nu1, params.nu2, params.lambda_2, params.omega_1, params.omega_2, params.phi, params.x_0, params.Lx, params.Ly, params.nuc_radius, params.theta_0];

%Choose tolerance sequence
%option_a = 1; %1 gives euclidean distance and mfpt etc. 0 gives spatial distribution and kl div etc.

%Generate fake data
%calculate appropriate summary statistic - we choose MFPT
%set prior
prior_params = [1.16, 0.8, 0.11, 0.42, 0.84, 0.58, 0.5, 0];
prior_sigma = [0.8, 0.8, 0.8]; %sd of gaussian or spread around mean of uniform
p_indices = [1, 4, 6];
par_params = prior_params;

abc_theta = zeros(N,length(p_indices));
fprintf('Generation 1 begins\n');

for i=1:N
    %initialise greater than tolerance
    %Uniform prior
    rr = rand(1,3);
    
    %simulate parameters from the prior
    par_params(p_indices) = prior_params(p_indices)+prior_sigma.*(rr-0.5);
    
    %simulate data using these model parameters
    %Calculate summary statistic (MFPT)
    abc_theta(i,:) = par_params(p_indices);
end

% figure(my_seed+1);
% subplot(3,1,1);
% plot(abc_theta(:,1),abc_theta(:,2),'o');
% hold all
% grid on
% plot(real_params(p_indices(1)),real_params(p_indices(2)),'rx','MarkerSize',12);
% set(gca, 'fontsize',14);
% xlabel('param1');
% ylabel('param2');
% subplot(3,1,2);
% plot(abc_theta(:,1),abc_theta(:,3),'o');
% hold all
% grid on
% plot(real_params(p_indices(1)),real_params(p_indices(3)),'rx','MarkerSize',12);
% set(gca, 'fontsize',14);
% xlabel('param1');
% ylabel('param3');
% subplot(3,1,3);
% plot(abc_theta(:,2),abc_theta(:,3),'o');
% hold all
% grid on
% plot(real_params(p_indices(2)),real_params(p_indices(3)),'rx','MarkerSize',12);
% set(gca, 'fontsize',14);
% xlabel('param2');
% ylabel('param3');

%entropy gives measure of difference from uniform distn
entropy = calculate_entropy(abc_theta,prior_params,prior_sigma,p_indices);

abc_weights = ones(size(abc_theta))/N;
%post-processing to check quality of posterior
[~, quality] = multi_dim_bin_posterior(abc_theta,abc_weights,prior_params,real_params,prior_sigma(1),p_indices,1);
%find 95% interval for posterior
box = zeros(length(p_indices),2);
contained_in_pred_interval = 1;
for k=1:3
    box(k,:) = [quantile(abc_theta(:,k),0.025),quantile(abc_theta(:,k),1-0.025)];
    contained_in_pred_interval = contained_in_pred_interval*(real_params(p_indices(k))>box(k,1))*(real_params(p_indices(k))<box(k,2));
end
box

% fname = sprintf('ABC_APMC_output%d.txt',my_seed);
% fileID = fopen(fname,'w');
% fprintf(fileID,'%f \n',abc_theta);
% fclose('all');

toc
end

