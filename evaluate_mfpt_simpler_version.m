function evaluate_mfpt_simpler_version(plot_option)
%last edit 15/6/15
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

param_vec = [0.1,0.5,1,5,10];%0.4:0.05:0.75; 0.5:0.02:0.8;
my_length = length(param_vec);
t_max = 200;
%my_length = 1;
mean_anchor_storage = zeros(my_length,1);
sd_anchor_storage = zeros(my_length,1);
mean_jumps_storage = zeros(my_length,1);
sd_jumps_storage = zeros(my_length,1);
mean_distances_storage = zeros(my_length,1);
sd_distances_storage = zeros(my_length,1);

num_particles = 100;
    

for k=1:my_length
    
    params.nu1 = param_vec(k);
    params.nu2 = param_vec(k)/2;
    %parfor? perhaps if it took longer might be needed
    anchored = zeros(num_particles,1);
    anchor_times = zeros(num_particles,1);
    num_jumps = zeros(num_particles,1);
    jump_distances = zeros(num_particles,1);
    
    for j=1:num_particles
        [anchored(j), anchor_times(j), ~, ~, pathx, pathy,~] = velocityjump2D_with_nucleus(t_max, params, 1, 2, 0);
        num_jumps(j) = length(pathx);
        jump_distances(j) = median(sqrt(diff(pathx).^2+diff(pathy).^2));
    end
    mean_anchor_time = mean(anchor_times);
    alt = sqrt(1/(num_particles-1)*sum((anchor_times-mean_anchor_time).^2));
    sd_anchor_time = std(anchor_times);
    
    mean_num_jumps = mean(num_jumps); %use mean number of jumps in a single passage as another
    sd_num_jumps = std(num_jumps);
    mean_jump_distances = mean(jump_distances); %use mean jump distance as final summary statistic
    sd_jump_distances = std(jump_distances);
    
    mean_anchor_storage(k) = mean_anchor_time;
    sd_anchor_storage(k) = sd_anchor_time;
    mean_jumps_storage(k) = mean_num_jumps;
    sd_jumps_storage(k) = sd_num_jumps;
    mean_distances_storage(k) = mean_jump_distances;
    sd_distances_storage(k) = sd_jump_distances;
end
mean_anchor_storage
sd_anchor_storage
if plot_option
    figure;
    subplot(3,1,1)
    errorbar(log10(param_vec),mean_anchor_storage,sd_anchor_storage,'linewidth',3)
    set(gca, 'fontsize',24);
    xlabel('log(\nu_1)');
    ylabel('MFPT');
    grid on
%    axis([param_vec(1)-0.05,param_vec(end)+0.05, 0, 3600*t_max]);

    %figure;
    subplot(3,1,2)
    errorbar(log10(param_vec),mean_jumps_storage,sd_jumps_storage,'linewidth',3)
    set(gca, 'fontsize',24);
    xlabel('log(\nu_1)');
    ylabel('Number of Jumps');
    grid on
    
    %figure;
    subplot(3,1,3)
    errorbar(log10(param_vec),mean_distances_storage,sd_distances_storage,'linewidth',3)
    set(gca, 'fontsize',24);
    xlabel('log(\nu_1)');
    ylabel('Jump Distance');
    grid on

end
toc(total)

end