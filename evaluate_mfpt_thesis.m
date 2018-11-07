function evaluate_mfpt_thesis()
%last edit 11/3/16
%created 6/5/15 harrison
%runs velocityjump3D_Nucleus which is a velocity jump process for a single particle
%runs this many times for many particles and returns simulated estimates of
%mfpt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath plotSpread/;
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

num_particles = 128;
my_length = 9;
param_vec = linspace(-1,1,my_length);
phi_vec = linspace(0.5,0.7,my_length); %0.5:0.02:0.8; [1.16/2,1.16];
param_indices = 1:5;
J = length(param_indices);
mfpt_storage = zeros(num_particles,my_length,J);
param_storage = zeros(6,my_length,J);


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
            params.lambda_2 = 10.^param_vec(k);            
        elseif j==3
            params.omega_1 = 10.^param_vec(k);
        elseif j==4
            params.omega_2 = 10.^param_vec(k);
        elseif j==5
            params.phi = phi_vec(k);  %0.5 + 0.5*exp(5*param_vec(k))/(1+exp(5*param_vec(k)));
        end
        anchored = zeros(num_particles,1);
        anchor_times = zeros(num_particles,1);
        num_jumps = zeros(num_particles,1);
        jump_distances = zeros(num_particles,1);
        param_storage(:,k,j) = [params.nu1, params.nu2, params.lambda_2, params.omega_1, params.omega_2, params.phi];
        [params.nu1, params.nu2, params.lambda_2, params.omega_1, params.omega_2, params.phi]
        
        parfor i=1:num_particles
            [anchored(i), anchor_times(i), ~, path, ~] = velocityjump3D_with_nucleus_and_new_BCs(params);
            num_jumps(i) = length(path(:,1));
            jump_distances(i) = median(sqrt(diff(path(:,1)).^2+diff(path(:,2)).^2+diff(path(:,3)).^2));
        end
%         mean_anchor_time = mean(anchor_times);
%         alt = sqrt(1/(num_particles-1)*sum((anchor_times-mean_anchor_time).^2));
%         sd_anchor_time = std(anchor_times);
%         
%         mean_num_jumps = mean(num_jumps); %use mean number of jumps in a single passage as another
%         sd_num_jumps = std(num_jumps);
%         mean_jump_distances = mean(jump_distances); %use mean jump distance as final summary statistic
%         sd_jump_distances = std(jump_distances);
        
        mfpt_storage(:,k,j) = anchor_times;
%         sd_anchor_storage(j,k) = sd_anchor_time;
%         mean_jumps_storage(j,k) = mean_num_jumps;
%         sd_jumps_storage(j,k) = sd_num_jumps;
%         mean_distances_storage(j,k) = mean_jump_distances;
%         sd_distances_storage(j,k) = sd_jump_distances;
    end
end
mfpt_storage = mfpt_storage/60^2; %scale as in hours for plotting
%sd_anchor_storage = sd_anchor_storage/60^2;
toc(total);
%output_matrix = [param_storage, mean_anchor_storage(:), sd_anchor_storage(:)];
%csvwrite('mfpt_sensitivity_to_R.csv',output_matrix);

if plot_option
    figure;
    plotSpread(mfpt_storage(:,:,1), ...
    'distributionIdx',reshape(repmat(param_vec,num_particles,1),[],1), ...
    'showMM',2, 'distributionColors','k');
    hold on;
    h = refline([0,0.3]);
    h.LineStyle = '--';
    h.Color='g';
    h.LineWidth=2;    
    set(gca, 'fontsize',24);
    xlabel('$\log_{10}(\nu) (\mu \textrm{m} \, \textrm{hr}^{-1})$','interpreter','latex');
    ylabel('MFPT \, (hr)','interpreter','latex');
    grid on;
    set(gca, 'XTick',[-1,-0.5,0,0.5,1]);
    set(gca, 'XTickLabel',[-1,-0.5,0,0.5,1]);
    fname = 'Figures_for_writeup/MFPT_nu_thesis';
    print(fname,'-depsc');
    %%%%%%%%%%%%%%%%%
    figure;
    plotSpread(mfpt_storage(:,:,2), ...
    'distributionIdx',reshape(repmat(param_vec,num_particles,1),[],1), ...
    'showMM',2, 'distributionColors','k');
    hold on;
    h = refline([0,0.3]);
    h.LineStyle = '--';
    h.Color='g';
    h.LineWidth=2;    
    set(gca, 'fontsize',24);
    xlabel('$\log_{10}(\lambda) (\textrm{hr}^{-1})$','interpreter','latex');
    ylabel('MFPT \, (hr)','interpreter','latex');
    grid on;
    set(gca, 'XTick',[-1,-0.5,0,0.5,1]);
    set(gca, 'XTickLabel',[-1,-0.5,0,0.5,1]);
    fname = 'Figures_for_writeup/MFPT_lambda_thesis';
    print(fname,'-depsc');
    %%%%%%%%%%%%%%
    figure;
    plotSpread(mfpt_storage(:,:,3), ...
    'distributionIdx',reshape(repmat(param_vec,num_particles,1),[],1), ...
    'showMM',2, 'distributionColors','k'); 
    hold on;
    h = refline([0,0.3]);
    h.LineStyle = '--';
    h.Color='g';
    h.LineWidth=2;
    set(gca, 'fontsize',24);
    xlabel('$\log_{10}(\omega_1) (\textrm{hr}^{-1})$','interpreter','latex');
    ylabel('MFPT \, (hr)','interpreter','latex');
    grid on;
    set(gca, 'XTick',[-1,-0.5,0,0.5,1]);
    set(gca, 'XTickLabel',[-1,-0.5,0,0.5,1]);
    fname = 'Figures_for_writeup/MFPT_omega_1_thesis';
    print(fname,'-depsc');
    %%%%%%%%%%%%%%%%%
    figure;
    plotSpread(mfpt_storage(:,:,4), ...
    'distributionIdx',reshape(repmat(param_vec,num_particles,1),[],1), ...
    'showMM',2, 'distributionColors','k'); 
    hold on;
    h = refline([0,0.3]);
    h.LineStyle = '--';
    h.Color='g';
    h.LineWidth=2;
    set(gca, 'fontsize',24);
    xlabel('$\log_{10}(\omega_2) (\textrm{hr}^{-1})$','interpreter','latex');
    ylabel('MFPT \, (hr)','interpreter','latex');
    grid on;
    set(gca, 'XTick',[-1,-0.5,0,0.5,1]);
    set(gca, 'XTickLabel',[-1,-0.5,0,0.5,1]);
    fname = 'Figures_for_writeup/MFPT_omega_2_thesis';
    print(fname,'-depsc');
    %%%%%%%%%%%%%%%%%
    figure;
    plotSpread(mfpt_storage(:,:,5), ...
    'distributionIdx',reshape(repmat(phi_vec,num_particles,1),[],1), ...
    'showMM',2, 'distributionColors','k');
    hold on;
    h = refline([0,0.3]);
    h.LineStyle = '--';
    h.Color='g';
    h.LineWidth=2;    
    set(gca, 'fontsize',24);
    xlabel('$\phi$','interpreter','latex');
    ylabel('MFPT \, (hr)','interpreter','latex');
    grid on;
    set(gca, 'XTick',min(phi_vec):0.1:max(phi_vec));
    set(gca, 'XTickLabel',min(phi_vec):0.1:max(phi_vec));    
    fname = 'Figures_for_writeup/MFPT_phi_thesis';
    print(fname,'-depsc');
end

end
