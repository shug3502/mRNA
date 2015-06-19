function [params, real_params, accepted_proportion, p_accept_min, ...
    q_estimate_fake, prior_params, p_indices, prior_sigma, par_params] = initialise_parameters(data_to_read,option_a)


%Fake parameters
    params.nu1 = 1.16; %speed of RNP complex under active transport [zimyanin et al 2008]
    params.nu2 = 0.80; %ratio between speed for active transport vs diffusion [zimyanin et al 2008]
    params.lambda_1=0;   %1/0.13; %transition rate =7.69 [zimyanin et al 2008]
    params.lambda_2 = 0.11;
    params.omega_1= 0.42;    %1/6*(num_modes>1); %rate of falling off the microtubule [zimyanin et al 2008] since average track length 2.4 - 2.8 microns -> average jump for 6s -> rate 1/6
    params.omega_2 = 0.84;
    params.phi = 0.58; %percentage of microtubules in posterior direction for biased angle distn [parton et al 2011]
    params.x_0=0.5;  %Initially in first compartment, ie. at NPC
    params.Lx = 52; %length of cell in x direction
    params.Ly = 37; %in y direction
    params.nuc_radius = 10; %radius of nucleus
    params.theta_0 = 0; %initial angle is 0

real_params = [params.nu1, params.nu2, params.lambda_2, params.omega_1, params.omega_2, params.phi, params.x_0, params.Lx, params.Ly, params.nuc_radius, params.theta_0];

%Choose tolerance sequence
accepted_proportion = 0.1; %alpha

p_accept_min = 0.05;
%option_a : 1 gives euclidean distance and mfpt etc. 0 gives spatial distribution and kl div etc.

if data_to_read
%read in data from files
fprintf('reading in data from file');
q_estimate_fake = read_in_particle_coords(9,[10,10]) %arguments are the file_number and the pixel sizes

else
%Generate fake data
%calculate appropriate summary statistic - we choose MFPT
fprintf('generating in silico data');
q_estimate_fake = summary_statistic_calculator(params,1000,0,option_a)
end
%create while loop

%set prior
prior_params = [5, 5, 5, 1.5, 1.5, 0.5, 0.5, 0];
%prior_params = real_params;
p_indices = 1:6;
%prior_params(p_indices) = [5,1.5,0.5];
prior_sigma = [10, 10, 10, 3, 3, 1]; 
% prior_params = [1.16, 0.8, 0.42, 0.42, 0.84, 0.58, 0.5, 0];
par_params = prior_params;

end