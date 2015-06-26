function q_averaged = run_multiple_VJ_with_production(repeats,with_plot)

    input_time = 25;
    with_anchoring = 1;
    num_modes = 1;
    %Now with data on nurse cells from Alex Davidson
    params.nu1 = 1.16; %speed of RNP complex under active transport [zimyanin et al 2008]
    params.nu2 = 0.80; %ratio between speed for active transport vs diffusion [zimyanin et al 2008]
    params.lambda_1=0;   %1/0.13; %transition rate =7.69 [zimyanin et al 2008]
    params.lambda_2 = 0.11;
    params.omega_1= 0.42;    %1/6*(num_modes>1); %rate of falling off the microtubule [zimyanin et al 2008] since average track length 2.4 - 2.8 microns -> average jump for 6s -> rate 1/6
    params.omega_2 = 0.84;
    params.phi = 0.58; %percentage of microtubules in posterior direction for biased angle distn [parton et al 2011]
    params.gamma = 0.08; %rate of prouction from the nucleus
    params.Lx = 52; %length of cell in x direction
    params.Ly = 37; %in y direction
    params.nuc_radius = 10; %radius of nucleus
    params.theta_0 = 0; %initial angle is 0
    
    delx=1;
    time_vec = 0:0.5:(input_time-1);
    l_t = length(time_vec);
    q_averaged = zeros(params.Lx/delx+1,l_t);
    
    for i=1:repeats
    [q_raw, ~] = velocityjump2D_with_production(input_time, params, with_anchoring, num_modes, 0);
    q_averaged = q_averaged + q_raw;
    end
    q_averaged = q_averaged./repeats;
if with_plot

    figure;
subplot(2,1,1)
imagesc(time_vec(1:9), 0:params.Lx, flipud(log(q_averaged(:,1:9)))) %plot initial period
colormap(gray)
set(gca, 'fontsize', 20);
yticklabels = 0:10:52;
yticks = linspace(1, size(q_averaged, 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)))
xlabel('Time')
ylabel('Position')
subplot(2,1,2)
imagesc(time_vec, 0:params.Lx, flipud(log(q_averaged))) % plot whole time course
colormap(gray)
set(gca, 'fontsize', 20);
yticklabels = 0:10:52;
yticks = linspace(1, size(q_averaged, 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)))
xlabel('Time')
ylabel('Position')

file_str = sprintf('Figures_for_writeup/Averaged_spatial_distn_production_gamma_%f',params.gamma);
print(file_str,'-deps');
print(file_str,'-dpng');
end

end