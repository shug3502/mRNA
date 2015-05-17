function evaluate_jump_lengths
%last edit 17/5/15
%created 17/5/15 harrison
%runs velocityjump2D_ModesInput which is a velocity jump process for a single particle
%runs this many times for many particles and returns simulated estimates of
%median jump lengths
total=tic;
close all

params.nu1 = 0.6;
params.nu2 = 0.1;
params.lambda = 1/0.13;
params.omega = 1/6;
params.phi = 0.58;
params.x_0 = 0.5;
params.theta_0 = 0;

t_max = 2;
phi_vec = 0.48:0.01:0.62;
my_length = length(phi_vec);
mean_num_jumps = zeros(my_length,1);
sd_storage = zeros(my_length,2);
mean_jump_distances = zeros(my_length,1);

for k=1:my_length
    
    num_particles = 20;
    %parfor? perhaps if it took longer might be needed
    num_jumps = zeros(num_particles,1);
    jump_distances = zeros(num_particles,1);
    
    params.phi = phi_vec(k);
    for j=1:num_particles
        [~, ~, ~, ~, pathx, ~] = velocityjump2D_ModesInput(t_max, params, 1, 0, 1, 0, 4*j);
        num_jumps(j) = length(pathx);
        jump_distances(j) = median(diff(pathx));
    end
    mean_num_jumps(k) = mean(num_jumps); %use mean number of jumps in a single passge as another
    mean_jump_distances(k) = mean(jump_distances); %use mean jump distance as final summary statistic
    sd_storage(k,:) = [std(num_jumps),std(jump_distances)];
end

sd_storage(:,1)
sd_storage(:,2)

    figure;
    subplot(2,1,1)
    errorbar(phi_vec,mean_num_jumps,sd_storage(:,1),'linewidth',3)
    set(gca, 'fontsize',16);
    xlabel('\phi');
    ylabel('mean # of jumps');
    grid on
    %axis([phi_vec(1)-0.5,phi_vec(end)+0.5, 0, 3600*t_max]);
    subplot(2,1,2)
    errorbar(phi_vec,mean_jump_distances,sd_storage(:,2),'r','linewidth',3)
    set(gca, 'fontsize',16);
    xlabel('\phi');
    ylabel('mean of median jump distances');
    grid on
    %axis([phi_vec(1)-0.5,phi_vec(end)+0.5, 0, 1]);

toc(total)

end