function [q_estimate, mfpt] = summary_statistic_calculator_production(par_params, useParamsVec)
if useParamsVec
    params.input_time = 1;
    params.with_anchoring = 1;
    params.num_modes = 1;
    params.with_plot = 0;
    params.nu1 = par_params(1); %speed of RNP complex under active transport [zimyanin et al 2008]
    params.nu2 = par_params(2); %ratio between speed for active transport vs diffusion [zimyanin et al 2008]
    params.lambda_1=0;   %1/0.13; %transition rate =7.69 [zimyanin et al 2008]
    params.lambda_2 = par_params(3);
    params.omega_1= par_params(4);    %1/6*(num_modes>1); %rate of falling off the microtubule [zimyanin et al 2008] since average track length 2.4 - 2.8 microns -> average$
    params.omega_2 = par_params(5);
    params.phi = par_params(6); %percentage of microtubules in posterior direction for biased angle distn [parton et al 2011]
    params.x_0= par_params(7);  %Initially in first compartment, ie. at NPC
    params.Lx = 52; %length of cell in x direction
    params.Ly = 37; %in y direction
    params.Lz = 37;
    params.nuc_radius = 10; %radius of nucleus
    params.rc_width = 1;
    params.ellipsoid_boundary = 1;
    params.gamma=0.1; %rate of production
else
    params = par_params;
end

t_end = params.input_time;
N=ceil(4*t_end*60^2*params.gamma);
%    delx = 1; %bin size;
L=params.Lx;
time_vec = 0:0.05:params.input_time;
t=time_vec*60^2;
l_t = length(time_vec);
jumps = zeros(N,10^3);
    q_estimate_spatial = zeros(3,l_t);
xpos = zeros(N,10^3);
%    xpos_discrete_time = zeros(N,10^3);
%    anchor_times = zeros(num_particles,1);
%    num_jumps = zeros(num_particles,1);
%    jump_speeds = zeros(num_particles,1);
%    lambda = zeros(num_particles,1);
%    mt_bias = zeros(num_particles,1);

time = 0;
num_particles=0;
while time < t_end
    num_particles=num_particles+1; %count num particles produced
    params.input_time = t_end-time/60^2;
    [is_anchored, ~,~, path, jump_times] = velocityjump3D_with_nucleus_and_new_BCs(params);
    xpos_temp = path(:,1);
    jump_temp = jump_times + time;
    jumps(num_particles,1:length(jump_temp)) = jump_temp;
    xpos(num_particles,1:length(xpos_temp)) = xpos_temp;
    for w=1:l_t
        find_index = find(t(w)<jumps(num_particles,:),1,'first');
        if isempty(xpos(num_particles,find_index)) && is_anchored
            xpos_discrete_time(num_particles,w) = L;
        else
            if find_index>1
                ddir = xpos(num_particles,find_index) - xpos(num_particles,find_index-1);
                xpos_discrete_time(num_particles,w) =  min((xpos(num_particles,find_index-1)+params.nu1*ddir/norm(ddir)*(t(w)-jump_times(find_index-1))),L);
            end
        end
    end
    
    %update time after one reaction (production event)
    tor = -1/params.gamma*log(rand(1));
    time = time + tor;
end
% error('stop now');
%
%     for j=1:N
%         [~, a_time, ~, xpos_temp, jump_temp] = velocityjump3D_with_nucleus_and_new_BCs(params);
% 	anchor_times(j) = a_time;
% 	path = xpos_temp;
% 	jump_times = jump_temp;
% 	num_jumps(j) = length(jump_times);
%         jump_speeds(j) = mean(sqrt(diff(path(:,1)).^2+diff(path(:,2)).^2+diff(path(:,3)).^2)./diff(jump_times)); %$
%         lambda(j) = 1/mean(diff(jump_times));
%         angle = acos(diff(path(:,1))./sum(diff(path).^2,2));
%         towards_posterior = (angle>pi/2).*(angle<3*pi/2);
%         mt_bias(j) = 1 - mean(towards_posterior);
%
%
% 	%NB jump_temp are the jump times and xpos_temp is the path
%         jumps(j,1:length(jump_temp)) = jump_temp;
%         xpos(j,1:size(xpos_temp,1)) = xpos_temp(:,1);
%         for w=1:l_t
%             if isempty(xpos(j,find(t(w)<jumps(j,:),1,'first')))
%                 xpos_discrete_time(j,w) = L/2;
%             else
%                 xpos_discrete_time(j,w) =  min(xpos(j,find(t(w)<jumps(j,:),1,'first')),L/2);
%             end
%         end
%     end
    for w=1:l_t
	%instead of splitting into bins, just store quantiles of data
        q_estimate_spatial(:,w) = quantile(xpos_discrete_time(:,w),3);
    end
%     mfpt = mean(anchor_times);  % use mean first passge time as one summary statistic
%     mean_num_jumps = mean(num_jumps); %use mean number of jumps in a single passage as another
%     mean_jump_speed = mean(jump_speeds); %use mean jump distance as final summary statistic
%     mean_lambda = mean(lambda);
%     mean_mt_bias = mean(mt_bias);
%     q_estimate_path = [mean_jump_speed; mean_lambda; mean_mt_bias];
 q_estimate = [num_particles; reshape(q_estimate_spatial,[],1)];
% q_estimate = [reshape(q_estimate_spatial,[],1); q_estimate_path];

end
