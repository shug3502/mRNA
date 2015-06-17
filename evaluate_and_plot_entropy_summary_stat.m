function evaluate_and_plot_entropy_summary_stat(plot_option)
%last edit 17/6/15
%created 16/5/15 harrison
%runs velocityjump2D_Nucleus which is a velocity jump process for a single particle
%runs this many times for many particles and returns simulated estimates of
%entropy
total=tic;
close all

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

param_vec = 0.4:0.05:0.8; %[0.001, 0.01,0.1,0.5,1,5,10]; 0.5:0.02:0.8;
my_length = length(param_vec);
%my_length = 1;

M=53;
q_storage = zeros(my_length,M);

num_particles = 100;
    

for k=1:my_length
    
    params.phi = param_vec(k);
    
    q_distn = summary_statistic_calculator(params,num_particles,0,0);
    q_storage(k,:) = q_distn;
end
if plot_option
    figure;
    imagesc((param_vec),0:52,q_storage')
    set(gca, 'fontsize',24);
    xlabel('\phi');
    ylabel('q');
    
end
toc(total)

end