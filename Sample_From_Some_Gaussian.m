function [abc_theta,abc_weights,entropy,quality,contained_in_pred_interval] = Sample_From_Some_Gaussian(N,option_a,my_seed)

%created 9/6/15
%last edit 9/6/15
%%sample a Gaussian approximation to a posterior about real parameters and
%%observe statistics
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
prior_sigma = [0.8, 0.8, 0.8]*2; %sd of gaussian or spread around mean of uniform
sigma = [0.4,0.4,0.4]/8;
p_indices = [1, 4, 6];
par_params = prior_params;

abc_theta = zeros(N,length(p_indices));
fprintf('Generation 1 begins\n');

for i=1:N
    %initialise greater than tolerance
    %Uniform prior
    rr = randn(1,3);
    
    %simulate parameters from the prior
    par_params(p_indices) = real_params(p_indices)+sigma.*rr;
    
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

abc_weights = ones(size(abc_theta))/N;
[entropy, quality, contained_in_pred_interval]= post_processing(abc_theta,abc_weights,prior_params,real_params,prior_sigma,p_indices);
toc
end

