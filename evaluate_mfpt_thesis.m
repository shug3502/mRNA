function evaluate_mfpt_thesis()
%last edit 11/3/16
%created 6/5/15 harrison
%runs velocityjump3D_Nucleus which is a velocity jump process for a single particle
%runs this many times for many particles and returns simulated estimates of
%mfpt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
total=tic;
close all
plot_option=1;
    params.input_time = 128;
    params.with_anchoring = 1;
    params.num_modes = 2;
    params.with_plot = 0;
    %Now with data on nurse cells from Alex Davidson
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
    params.Lz = 37; %in z direction
    params.nuc_radius = 10; %radius of nucleus
    params.theta_0 = 0; %initial angle is 0
    params.rc_width = 1;
    params.ellipsoid_boundary = 1;

param_vec = linspace(-1,1,10); %[0.1,0.5,1,5,10]; %0.4:0.05:0.75; 0.5:0.02:0.8; [1.16/2,1.16];
param_indices = [1,2,6];
my_length = length(param_vec);
J = length(param_indices);
mean_anchor_storage = zeros(J,my_length);
sd_anchor_storage = zeros(J,my_length);
mean_jumps_storage = zeros(J,my_length);
sd_jumps_storage = zeros(J,my_length);
mean_distances_storage = zeros(J,my_length);
sd_distances_storage = zeros(J,my_length);
param_storage = zeros(J*my_length,6);

num_particles = 100;
    
for j=1:J
    %reset changes to parameters
    params.nu1 = 1.16; %speed of RNP complex under active transport [zimyanin et al 2008]
    params.nu2 = 0.80; %ratio between speed for active transport vs diffusion [zimyanin et al 2008]
    params.lambda_2 = 0.11;
    params.omega_1= 0.42;    %1/6*(num_modes>1); %rate of falling off the microtubule [zimyanin et al 2008] since average track length 2.4 - 2.8 microns -> average jump for 6s -> rate 1/6
    params.omega_2 = 0.84;
    params.phi = 0.58;
    for k=1:my_length
        k
        if j==1 %may consider other parameters or combinations of parameters
            params.nu1 = 10.^param_vec(k);
            params.nu2 = (10.^param_vec(k))/2;
        elseif j==2
            params.omega_1 = 10.^param_vec(k);
        elseif j==3
            params.phi = 0.5 + 0.5*exp(5*param_vec(k))/(1+exp(5*param_vec(k)));
        end
    anchored = zeros(num_particles,1);
    anchor_times = zeros(num_particles,1);
    num_jumps = zeros(num_particles,1);
    jump_distances = zeros(num_particles,1);
    param_storage((j-1)*my_length + k,:) = [params.nu1, params.nu2, params.lambda_2, params.omega_1, params.omega_2, params.phi];
    [params.nu1, params.nu2, params.lambda_2, params.omega_1, params.omega_2, params.phi]
    
    parfor i=1:num_particles
	[anchored(i), anchor_times(i), ~, path, ~] = velocityjump3D_with_nucleus_and_new_BCs(params);
        num_jumps(i) = length(path(:,1));
        jump_distances(i) = median(sqrt(diff(path(:,1)).^2+diff(path(:,2)).^2+diff(path(:,3)).^2));
    end
    mean_anchor_time = mean(anchor_times);
    alt = sqrt(1/(num_particles-1)*sum((anchor_times-mean_anchor_time).^2));
    sd_anchor_time = std(anchor_times);
    
    mean_num_jumps = mean(num_jumps); %use mean number of jumps in a single passage as another
    sd_num_jumps = std(num_jumps);
    mean_jump_distances = mean(jump_distances); %use mean jump distance as final summary statistic
    sd_jump_distances = std(jump_distances);
    
    mean_anchor_storage(j,k) = mean_anchor_time;
    sd_anchor_storage(j,k) = sd_anchor_time;
    mean_jumps_storage(j,k) = mean_num_jumps;
    sd_jumps_storage(j,k) = sd_num_jumps;
    mean_distances_storage(j,k) = mean_jump_distances;
    sd_distances_storage(j,k) = sd_jump_distances;
end
end
mean_anchor_storage = mean_anchor_storage/60^2; %scale as in hours for plotting
sd_anchor_storage = sd_anchor_storage/60^2; 
toc(total);
%output_matrix = [param_storage, mean_anchor_storage(:), sd_anchor_storage(:)];
%csvwrite('mfpt_sensitivity_to_R.csv',output_matrix);

if plot_option
    figure;
    errorbar(10.^param_vec,mean_anchor_storage(1,:),sd_anchor_storage(1,:),'linewidth',3)
    hold on;
    h = refline([0,0.3]);
    h.LineStyle = '--';
    h.Color='g';
    set(gca, 'fontsize',24);
    xlabel('\nu (\mu m hr^{-1})');
    ylabel('MFPT');
    grid on
    %axis([0,12, 0,5000]);
    fname = 'Figures_for_writeup/MFPT_nu_thesis';
    print(fname,'-depsc');
%%%%%%%%%%%%%%%%%
    figure;
    errorbar(10.^param_vec,mean_anchor_storage(2,:),sd_anchor_storage(2,:),'linewidth',3);
    hold on;
    h = refline([0,0.3]);
    h.LineStyle = '--';
    h.Color='g';
    set(gca, 'fontsize',24);
    xlabel('\omega_1 (hr^{-1})');
    ylabel('MFPT');
    grid on
    %axis([0,12, 0,5000]);
    fname = 'Figures_for_writeup/MFPT_omega_1_thesis';
    print(fname,'-depsc');
%%%%%%%%%%%%%%%%%
    figure;
    errorbar(0.5 + 0.5*exp(5*param_vec)./(1+exp(5*param_vec)),mean_anchor_storage(3,:) ,sd_anchor_storage(3,:),'linewidth',3);
    hold on;
    h = refline([0,0.3]);
    h.LineStyle = '--';
    h.Color='g';
    set(gca, 'fontsize',24);
    xlabel('\phi');
    ylabel('MFPT');
    grid on
    %axis([0,12, 0,5000]);
    fname = 'Figures_for_writeup/MFPT_phi_thesis';
    print(fname,'-depsc');
end

end
