function evaluate_mfpt_simpler_version(plot_option)
%last edit 9/6/15
%created 6/5/15 harrison
%runs velocityjump2D_Nucleus which is a velocity jump process for a single particle
%runs this many times for many particles and returns simulated estimates of
%mfpt
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

t_max = 100;
my_length = 1;
mean_anchor_storage = zeros(my_length,1);
sd_anchor_storage = zeros(my_length,1);

num_particles = 500;

for k=1:my_length
    
    
    %parfor? perhaps if it took longer might be needed
    anchored = zeros(num_particles,1);
    anchor_times = zeros(num_particles,1);
    final_positions = zeros(num_particles,2);
    
    for j=1:num_particles
        [anchored(j), anchor_times(j), final_positions(j,1), final_positions(j,2), ~,~,~] = velocityjump2D_with_nucleus(t_max, params, 1,2, 0);
    end
    mean_anchor_time = mean(anchor_times);
    alt = sqrt(1/(num_particles-1)*sum((anchor_times-mean_anchor_time).^2));
    sd_anchor_time = std(anchor_times);
    
    mean_anchor_storage(k) = mean_anchor_time;
    sd_anchor_storage(k) = sd_anchor_time;
end
mean_anchor_storage
sd_anchor_storage
if plot_option
    x_0_vec=0.5;
    figure;
    errorbar(x_0_vec,mean_anchor_storage,sd_anchor_storage,'linewidth',3)
    set(gca, 'fontsize',16);
    xlabel('x_0');
    ylabel('MFPT');
    grid on
    axis([x_0_vec(1)-0.5,x_0_vec(end)+0.5, 0, 3600*t_max]);
end
toc(total)

end