function [q_estimate] = summary_statistic_calculator(par_params,num_particles,is_parallel,option_a_summary_statistic)

if option_a_summary_statistic %choose which type of summary statistic to use
    
    %runs velocityjump2D_ModesInput which is a velocity jump process for a single particle
    %runs this many times for many particles and returns simulated estimates of
    %mfpt
    
    if is_parallel
        params.nu1 = par_params(1); %speed of RNP complex under active transport [zimyanin et al 2008]
        params.nu2 = par_params(2); %ratio between speed for active transport vs diffusion [zimyanin et al 2008]
        params.lambda_1=0;   %1/0.13; %transition rate =7.69 [zimyanin et al 2008]
        params.lambda_2 = par_params(3);
        params.omega_1= par_params(4);    %1/6*(num_modes>1); %rate of falling off the microtubule [zimyanin et al 2008] since average track length 2.4 - 2.8 microns -> average jump for 6s -> rate 1/6
        params.omega_2 = par_params(5);
        params.phi = par_params(6); %percentage of microtubules in posterior direction for biased angle distn [parton et al 2011]
        params.x_0= par_params(7);  %Initially in first compartment, ie. at NPC
        params.Lx = 52; %length of cell in x direction
        params.Ly = 37; %in y direction
        params.nuc_radius = 10; %radius of nucleus
        params.theta_0 = par_params(8); %initial angle is 0
        
    else
        params = par_params;
    end
    non_infinite = 0;
    t_max = 100;  %if this is too small, some runs will not reach anchoring giving infinite mfpt _> rejection.
    % these are probably not an issue as they would give large mfpt anyway,
    % which would be rejected. But can increase this.
    
    anchor_times = zeros(num_particles,1);
    num_jumps = zeros(num_particles,1);
    jump_distances = zeros(num_particles,1);
    parfor j=1:num_particles
        [~, anchor_times(j), ~, ~, pathx, pathy, ~] = velocityjump2D_with_nucleus(t_max, params, 1, 2, 0);
        num_jumps(j) = length(pathx);
        jump_distances(j) = median(sqrt(diff(pathx).^2+diff(pathy).^2)); %median(abs(diff(pathx)))
    end
    mean_fp_time = mean(anchor_times);  % use mean first passge time as one summary statistic
    mean_num_jumps = mean(num_jumps); %use mean number of jumps in a single passage as another
    mean_jump_distances = mean(jump_distances); %use mean jump distance as final summary statistic
    q_estimate = [mean_fp_time; mean_num_jumps; mean_jump_distances];
    %non_infinite = 0;
    % if non_infinite
    % while isinf(mean_fp_time)
    %     t_max = 2*t_max;
    %     parfor j=1:num_particles
    %     [~, anchor_times(j), ~, ~, pathx, ~] = velocityjump2D_ModesInput(t_max, params, 1, 0, 1, 0, r_random*j);
    %     num_jumps(j) = length(pathx);
    %     end
    %     mean_fp_time = mean(anchor_times);
    %     mean_num_jumps = mean(num_jumps);
    % end
    % end
else %else use the spatial distribution and kl divergence
    
    %runs velocityjump2D_ModesInput which is a velocity jump process for a single particle
    %runs this many times for many particles and returns simulated estimates of
    %mfpt
    
    delx = 1; %bin size;
    if is_parallel
        params.nu1 = par_params(1); %speed of RNP complex under active transport [zimyanin et al 2008]
        params.nu2 = par_params(2); %ratio between speed for active transport vs diffusion [zimyanin et al 2008]
        params.lambda_1=0;   %1/0.13; %transition rate =7.69 [zimyanin et al 2008]
        params.lambda_2 = par_params(3);
        params.omega_1= par_params(4);    %1/6*(num_modes>1); %rate of falling off the microtubule [zimyanin et al 2008] since average track length 2.4 - 2.8 microns -> average jump for 6s -> rate 1/6
        params.omega_2 = par_params(5);
        params.phi = par_params(6); %percentage of microtubules in posterior direction for biased angle distn [parton et al 2011]
        params.x_0= par_params(7);  %Initially in first compartment, ie. at NPC
        params.Lx = 52; %length of cell in x direction
        params.Ly = 37; %in y direction
        params.nuc_radius = 10; %radius of nucleus
        params.theta_0 = par_params(8); %initial angle is 0
        
    else
        params = par_params;
    end
    
    N=num_particles;
    L=params.Lx;
    time_vec = (0.15:0.05:0.15)*0.5;
    t=time_vec*60^2;
    l_t = length(time_vec);
    jumps = zeros(N,10^3);
    q_estimate = zeros(L/delx+1,l_t);
    xpos = zeros(N,10^3);
    xpos_discrete_time = zeros(N,10^3);
    
    for j=1:N
        [~, ~, ~, ~, xpos_temp, ~, jump_temp] = velocityjump2D_with_nucleus(max(time_vec), params, 1, 2, 0);
        jumps(j,1:length(jump_temp)) = jump_temp;
        xpos(j,1:length(xpos_temp)) = xpos_temp;
        for w=1:l_t
            if isempty(xpos(j,find(t(w)<jumps(j,:),1,'first')))
                xpos_discrete_time(j,w) = L;
            else
                xpos_discrete_time(j,w) =  min(xpos(j,find(t(w)<jumps(j,:),1,'first')),L);
            end
        end
    end
    for w=1:l_t
        %split into bins
        [Num_in_bins,~] = histc(xpos_discrete_time(:,w),0:delx:L);
        if sum(Num_in_bins) ~= N
            xpos_discrete_time(:,w)
            ename = sprintf('ABC_errors.txt');
            fileIDerror = fopen(ename,'w');
            fprintf(fileIDerror,'%f \n',xpos_discrete_time);
            fclose('all');
            error('outside of 0:L');
        end
        q_estimate(:,w) = Num_in_bins/N/delx; %estimate of q at time T
    end
end

end