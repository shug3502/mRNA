function calculate_population_statistics(plot_option)
%last edit 8/5/15
%created 21/4/15 harrison
%runs velocityjump2D which is a velocity jump process for a single particle
%runs this many times for many particles and returns population level statistics
%edited to take parameters as a struct
total=tic;
close all

params.nu1 = 0.4;
params.nu2 = 0.08;
params.lambda = 1/0.13;
params.omega = 1/6;
params.x_0=0.5;
params.theta_0 = 0;

phi_vec = 48:65;
mean_anchor_storage = zeros(length(phi_vec),1);
sd_anchor_storage = zeros(length(phi_vec),1);
proportion_anchor_storage = zeros(length(phi_vec),1);
for k=1:length(phi_vec)
    
    params.phi = phi_vec(k)/100;
    num_particles = 10;
    %parfor? perhaps if it took longer might be needed
    anchored = zeros(num_particles,1);
    anchor_times = zeros(num_particles,1);
    final_positions = zeros(num_particles,2);
    for j=1:num_particles
        [anchored(j),anchor_times(j),final_positions(j,1),final_positions(j,2)]  = velocityjump2D_ModesInput(6, params, 1, 0, 2, plot_option, j);
    end
    proportion_anchored = sum(anchored)/length(anchored)
    mean_anchor_time = mean(anchor_times)
    sd_anchor_time = std(anchor_times)
    
    if plot_option
        figure;
        subplot(2,1,1)
        hist(final_positions(:,1))
        set(gca,'fontsize',16);
        xlabel('x position')
        ylabel('frequency')
        subplot(2,1,2)
        hist(final_positions(:,2))
        set(gca,'fontsize',16);
        xlabel('y position')
        ylabel('frequency')
    end
    
    mean_anchor_storage(k) = mean_anchor_time;
    sd_anchor_storage(k) = sd_anchor_time;
    proportion_anchor_storage(k) = proportion_anchored;
end

figure;
subplot(2,1,1)
errorbar(phi_vec,mean_anchor_storage,sd_anchor_storage,'linewidth',3)
set(gca, 'fontsize',16);
xlabel('\phi');
ylabel('MFPT');
grid on
axis([phi_vec(1),phi_vec(end), 0, 5000]);
subplot(2,1,2)
plot(phi_vec,proportion_anchor_storage,'r','linewidth',3)
set(gca, 'fontsize',16);
xlabel('\phi');
ylabel('Proportion anchored after 6hrs');
grid on
axis([phi_vec(1),phi_vec(end), 0, 1]);

toc(total)

end