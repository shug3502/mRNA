function evaluate_mfpt(plot_option,vary_x,vary_lambda)
%last edit 6/5/15
%created 6/5/15 harrison
%runs velocityjump2D_ModesInput which is a velocity jump process for a single particle
%runs this many times for many particles and returns simulated estimates of
%mfpt
total=tic;
close all

if vary_x*vary_lambda==1
    error('Cannot vary both lambda and x');
end
params.nu1 = 0.4;
params.nu2 = 0.1;
params.lambda = 1/0.13;
params.omega = 1/6;
params.phi = 0.58;

t_max = 1;
x_0 = 0.5;
x_0_vec = 0.5:4:20.5;
theta_0 = 0;
theta_0_vec = linspace(0,2*pi,7);
lambda_vec = 0.5:0.5:30;
my_length = (length(x_0_vec)*vary_x + length(theta_0_vec)*(1-vary_x))*(1-vary_lambda) + vary_lambda*length(lambda_vec);
mean_anchor_storage = zeros(my_length,1);
sd_anchor_storage = zeros(my_length,1);
proportion_anchor_storage = zeros(my_length,1);

for k=1:my_length
    
    num_particles = 20;
    %parfor? perhaps if it took longer might be needed
    anchored = zeros(num_particles,1);
    anchor_times = zeros(num_particles,1);
    final_positions = zeros(num_particles,2);
    if vary_x  %we vary the initial value of x_0
        params.x_0 = x_0_vec(k);
        params.theta_0 = theta_0;
        for j=1:num_particles
            [anchored(j), anchor_times(j), final_positions(j,1),final_positions(j,2)] = velocityjump2D_ModesInput(t_max, params, 1, 0, 1, plot_option, 4*j);
        end
    else %we vary the initial angle theta_0
        if vary_lambda
        params.x_0 = x_0;
        params.theta_0 = theta_0;
        params.lambda = lambda_vec(k);
        for j=1:num_particles
            [anchored(j), anchor_times(j), final_positions(j,1),final_positions(j,2)] = velocityjump2D_ModesInput(t_max, params, 1, 0, 1, plot_option, 4*j);
        end    
        else
        params.x_0 = x_0;
        params.theta_0 = theta_0_vec(k);
        for j=1:num_particles
            [anchored(j), anchor_times(j), final_positions(j,1),final_positions(j,2)] = velocityjump2D_ModesInput(t_max, params, 1, 0, 1, plot_option, 4*j);
        end
        end
    end
    proportion_anchored = sum(anchored)/length(anchored)
    mean_anchor_time = mean(anchor_times)
    sd_anchor_time = std(anchor_times)
    
    mean_anchor_storage(k) = mean_anchor_time;
    sd_anchor_storage(k) = sd_anchor_time;
    proportion_anchor_storage(k) = proportion_anchored;
end

if vary_x
    figure;
    subplot(2,1,1)
    errorbar(x_0_vec,mean_anchor_storage,sd_anchor_storage,'linewidth',3)
    set(gca, 'fontsize',16);
    xlabel('x_0');
    ylabel('MFPT');
    grid on
    axis([x_0_vec(1)-0.5,x_0_vec(end)+0.5, 0, 3600*t_max]);
    subplot(2,1,2)
    plot(x_0_vec,proportion_anchor_storage,'r','linewidth',3)
    set(gca, 'fontsize',16);
    xlabel('x_0');
    ylabel('Proportion anchored after 0.5hrs');
    grid on
    axis([x_0_vec(1)-0.5,x_0_vec(end)+0.5, 0, 1]);
else
    if vary_lambda
    figure;
    subplot(2,1,1)
    errorbar(lambda_vec,mean_anchor_storage,sd_anchor_storage,'linewidth',3)
    set(gca, 'fontsize',16);
    xlabel('\lambda');
    ylabel('MFPT');
    grid on
    axis([lambda_vec(1),lambda_vec(end), 0, 1000*t_max]);
    subplot(2,1,2)
    plot(lambda_vec,proportion_anchor_storage,'r','linewidth',3)
    set(gca, 'fontsize',16);
    xlabel('\lambda');
    ylabel('Proportion anchored after 0.5hrs');
    grid on
    axis([lambda_vec(1),lambda_vec(end), 0, 1]);
    else
    figure;
    subplot(2,1,1)
    errorbar(theta_0_vec,mean_anchor_storage,sd_anchor_storage,'linewidth',3)
    set(gca, 'fontsize',16);
    xlabel('\theta_0');
    ylabel('MFPT');
    grid on
    axis([theta_0(1),theta_0_vec(end), 0, 3600*t_max]);
    subplot(2,1,2)
    plot(theta_0_vec,proportion_anchor_storage,'r','linewidth',3)
    set(gca, 'fontsize',16);
    xlabel('\theta_0');
    ylabel('Proportion anchored after 0.5hrs');
    grid on
    axis([theta_0_vec(1),theta_0_vec(end), 0, 1]);
    end
end

toc(total)

end