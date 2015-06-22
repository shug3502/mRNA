function [abc_theta,abc_weights,entropy,quality,contained_in_pred_interval] = ABC_RS_knn_production(N,my_seed,data_to_read)

%created 22/6/15
%last edit 22/6/15
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
%Now applied to new model with production

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
rng(my_seed);
close all


%Fake parameters
    params.nu1 = 1.16; %speed of RNP complex under active transport [zimyanin et al 2008]
    params.nu2 = 0.80; %ratio between speed for active transport vs diffusion [zimyanin et al 2008]
    params.lambda_1=0;   %1/0.13; %transition rate =7.69 [zimyanin et al 2008]
    params.lambda_2 = 0.11;
    params.omega_1= 0.42;    %1/6*(num_modes>1); %rate of falling off the microtubule [zimyanin et al 2008] since average track length 2.4 - 2.8 microns -> average jump for 6s -> rate 1/6
    params.omega_2 = 0.84;
    params.phi = 0.58; %percentage of microtubules in posterior direction for biased angle distn [parton et al 2011]
    params.gamma = 0.01; %production
    params.Lx = 52; %length of cell in x direction
    params.Ly = 37; %in y direction
    params.nuc_radius = 10; %radius of nucleus
    params.theta_0 = 0; %initial angle is 0

real_params = [params.nu1, params.nu2, params.lambda_2, params.omega_1, params.omega_2, params.phi, params.gamma, params.Lx, params.Ly, params.nuc_radius, params.theta_0];

%Choose tolerance sequence
accepted_proportion = 0.1; %alpha

if data_to_read
%read in data from files
fprintf('reading in data from file');
q_estimate_fake = read_in_particle_coords(9,[10,10]) %arguments are the file_number and the pixel sizes

else
%Generate fake data
%calculate appropriate summary statistic - we choose MFPT
fprintf('generating in silico data');
q_estimate_fake = generate_summary_stat(5, real_params, 10)
end
%create while loop

%set prior
prior_params = [5, 5, 5, 1.5, 1.5, 0.5, 0.5, 0];
%prior_params = real_params;
p_indices = 1:7;
%prior_params(p_indices) = [5,1.5,0.5];
prior_sigma = [10, 10, 10, 3, 3, 1, 1]; 
% prior_params = [1.16, 0.8, 0.42, 0.42, 0.84, 0.58, 0.01, 0];
par_params = prior_params;




abc_theta = zeros(N,length(p_indices));
abc_weights = zeros(N,1);
abc_dist = zeros(N,1);
fprintf('Generation 1 begins\n');

for i=1:N
    if mod(i,10)==0
        i
    end
    %initialise greater than tolerance
    %Uniform prior
    rr = rand(1,length(p_indices));
    
    %simulate parameters from the prior
    par_params(p_indices) = prior_params(p_indices)+prior_sigma(p_indices).*(rr-0.5);
    
    %simulate data using these model parameters
    %Calculate summary statistic (MFPT)
    q_estimate_candidate = generate_summary_stat(5, par_params, 1);
    abc_dist(i) = distance_metric(q_estimate_candidate,q_estimate_fake,params.Lx,0); %distance of proposed S(x) from S(x_obs)
    %end    %repeat until N acceptances have been made
    
    abc_theta(i,:) = par_params(p_indices);
    abc_weights(i) = 1/N;
end

%keep now parameters with distances less than the accepted quantile
to_keep = (abc_dist <= quantile(abc_dist,accepted_proportion));
abc_theta = abc_theta(to_keep,:);
abc_weights = abc_weights(to_keep)./sum(abc_weights(to_keep));
abc_dist = abc_dist(to_keep);

if length(p_indices)>=3
figure(my_seed+1);
subplot(3,1,1);
plot(abc_theta(:,1),abc_theta(:,2),'o');
hold all
grid on
plot(real_params(p_indices(1)),real_params(p_indices(2)),'rx','MarkerSize',12);
set(gca, 'fontsize',14);
xlabel('param1');
ylabel('param2');
subplot(3,1,2);
plot(abc_theta(:,1),abc_theta(:,3),'o');
hold all
grid on
plot(real_params(p_indices(1)),real_params(p_indices(3)),'rx','MarkerSize',12);
set(gca, 'fontsize',14);
xlabel('param1');
ylabel('param3');
subplot(3,1,3);
plot(abc_theta(:,2),abc_theta(:,3),'o');
hold all
grid on
plot(real_params(p_indices(2)),real_params(p_indices(3)),'rx','MarkerSize',12);
set(gca, 'fontsize',14);
xlabel('param2');
ylabel('param3');
end

[entropy, quality, contained_in_pred_interval]= post_processing(abc_theta,abc_weights,prior_params,real_params,prior_sigma,p_indices);

toc
end

function q_summary = generate_summary_stat(input_time, vparams, repeats)
    params.nu1 = vparams(1); %speed of RNP complex under active transport [zimyanin et al 2008]
    params.nu2 = vparams(2); %ratio between speed for active transport vs diffusion [zimyanin et al 2008]
    params.lambda_1=0;   %1/0.13; %transition rate =7.69 [zimyanin et al 2008]
    params.lambda_2 = vparams(3);
    params.omega_1= vparams(4);    %1/6*(num_modes>1); %rate of falling off the microtubule [zimyanin et al 2008] since average track length 2.4 - 2.8 microns -> average jump for 6s -> rate 1/6
    params.omega_2 = vparams(5);
    params.phi = vparams(6); %percentage of microtubules in posterior direction for biased angle distn [parton et al 2011]
    params.gamma = vparams(7); %production
    params.Lx = 52; %length of cell in x direction
    params.Ly = 37; %in y direction
    params.nuc_radius = 10; %radius of nucleus
    params.theta_0 = 0; %initial angle is 0



delx=1;
q_summary = zeros(params.Lx/delx+1,(input_time-1)/0.5+1);
for i=1:repeats
    q_estimate_fake = velocityjump2D_with_production(input_time, params, 1, 2, 0);  %input time of 5 hrs
    q_summary = q_summary+q_estimate_fake;
end
q_summary=q_summary/repeats;
end

