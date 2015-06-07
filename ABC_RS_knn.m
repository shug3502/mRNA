function [abc_theta,abc_weights,entropy,quality,contained_in_pred_interval] = ABC_RS_knn(N,option_a,my_seed)

%created 4/6/15
%last edit 28/5/15
%AS simple as ABC
%Based on Simplest_ABC_with_moments
%added in knn selection of theta for parameters rather than less than a
%specific delta value
%implements a version of ABC via population mc to fit parameters to
%some 'data' which is first simulated from the model
%Model used is velocity jump process with 1 mode
%see Turner 2012 for ABC population monte carlo tutorial
%Should be extendable to adaptive popMc-ABC
%See Lenormand 2013 for APMC
%Should be edited to now use new code for vel jump with nucleus

%Simple rejection sampling ABC on complicated model to match other code for
%APMC etc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
rng(my_seed);
close all

%Fake parameters
params.nu1 = 1.0; %speed of RNP complex under active transport [zimyanin et al 2008]
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

real_params = [params.nu1, params.nu2, params.lambda_2, params.omega_1, params.omega_2, params.phi, params.x_0, params.Lx, params.Ly, params.nuc_radius, params.theta_0];

%Choose tolerance sequence
accepted_proportion = 0.5; %0.5 %alpha
%At t=1 for first generation
%N=500;

%option_a = 1; %1 gives euclidean distance and mfpt etc. 0 gives spatial distribution and kl div etc.

%Generate fake data
%calculate appropriate summary statistic - we choose MFPT
q_estimate_fake = summary_statistic_calculator(params,1000,0,option_a)

%create while loop

%set prior
prior_params = [1.16, 0.8, 0.11, 0.42, 0.84, 0.58, 0.5, 0];
prior_sigma = [0.8, 0.8, 0.8]; %sd of gaussian or spread around mean of uniform
p_indices = [1, 4, 6];
par_params = prior_params;

abc_theta = zeros(N,length(p_indices));
abc_weights = zeros(N,1);
abc_dist = zeros(N,1);
fprintf('Generation 1 begins\n');

for i=1:N
    %initialise greater than tolerance
    %Uniform prior
    rr = rand(1,3);
    
    %simulate parameters from the prior
    par_params(p_indices) = prior_params(p_indices)+prior_sigma.*(rr-0.5);
    
    %simulate data using these model parameters
    %Calculate summary statistic (MFPT)
    q_estimate_candidate = summary_statistic_calculator(par_params,20,1,option_a);
    abc_dist(i) = distance_metric(q_estimate_candidate,q_estimate_fake,params.Lx,option_a); %distance of proposed S(x) from S(x_obs)
    %end    %repeat until N acceptances have been made
    
    abc_theta(i,:) = par_params(p_indices);
    abc_weights(i) = 1/N;
end

%keep now parameters with distances less than the accepted quantile
to_keep = (abc_dist <= quantile(abc_dist,accepted_proportion));
abc_theta = abc_theta(to_keep,:);
abc_weights = abc_weights(to_keep)./sum(abc_weights(to_keep));
abc_dist = abc_dist(to_keep);

figure(my_seed+1);
subplot(3,1,1);
plot(abc_theta(:,1),abc_theta(:,2),'o');
hold all
grid on
plot(real_params(p_indices(1)),real_params(p_indices(2)),'rx','MarkerSize',12);
set(gca, 'fontsize',14);
xlabel('param1');
ylabel('param2');
subplot(3,1,2);
plot(abc_theta(:,1),abc_theta(:,3),'o');
hold all
grid on
plot(real_params(p_indices(1)),real_params(p_indices(3)),'rx','MarkerSize',12);
set(gca, 'fontsize',14);
xlabel('param1');
ylabel('param3');
subplot(3,1,3);
plot(abc_theta(:,2),abc_theta(:,3),'o');
hold all
grid on
plot(real_params(p_indices(2)),real_params(p_indices(3)),'rx','MarkerSize',12);
set(gca, 'fontsize',14);
xlabel('param2');
ylabel('param3');

%entropy gives measure of difference from uniform distn
entropy = calculate_entropy(abc_theta,prior_params,prior_sigma,p_indices);
abc_theta

%post-processing to check quality of posterior
[~, quality] = multi_dim_bin_posterior(abc_theta,abc_weights,prior_params,real_params,prior_sigma(1),p_indices,1);
%find 95% interval for posterior
box = zeros(length(p_indices),2);
contained_in_pred_interval = 1;
for k=1:3
    box(k,:) = [quantile(abc_theta(:,k),0.025),quantile(abc_theta(:,k),1-0.025)];
    contained_in_pred_interval = contained_in_pred_interval*(real_params(p_indices(k))>box(k,1))*(real_params(p_indices(k))<box(k,2));
end
box

% fname = sprintf('ABC_APMC_output%d.txt',my_seed);
% fileID = fopen(fname,'w');
% fprintf(fileID,'%f \n',abc_theta);
% fclose('all');

toc
end

%use mfpt, num_jumps and jump distances as summary statistics
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
    %non_infinite = 0;
    t_max = 1;  %if this is too small, some runs will not reach anchoring giving infinite mfpt _> rejection.
    % these are probably not an issue as they would give large mfpt anyway,
    % which would be rejected. But can increase this.
    
    anchor_times = zeros(num_particles,1);
    num_jumps = zeros(num_particles,1);
    jump_distances = zeros(num_particles,1);
    for j=1:num_particles
        [~, anchor_times(j), ~, ~, pathx, pathy, ~] = velocityjump2D_with_nucleus(t_max, params, 1, 2, 0);
        num_jumps(j) = length(pathx);
        jump_distances(j) = median(sqrt(diff(pathx).^2+diff(pathy).^2)); %median(abs(diff(pathx)))
    end
    mean_fp_time = mean(anchor_times);  % use mean first passge time as one summary statistic
    mean_num_jumps = mean(num_jumps); %use mean number of jumps in a single passage as another
    mean_jump_distances = mean(jump_distances); %use mean jump distance as final summary statistic
    q_estimate = [mean_fp_time; mean_num_jumps; mean_jump_distances];
%     non_infinite = 1;
%     if non_infinite
%         while isinf(mean_fp_time)
%             t_max = 2*t_max
%             for j=1:num_particles
%                 [~, anchor_times(j), ~, ~, pathx, pathy, ~] = velocityjump2D_with_nucleus(t_max, params, 1, 2, 0);
%                 num_jumps(j) = length(pathx);
%                 jump_distances(j) = median(sqrt(diff(pathx).^2+diff(pathy).^2)); %median(abs(diff(pathx)))
%             end
%             mean_fp_time = mean(anchor_times);
%             mean_num_jumps = mean(num_jumps);
%             mean_jump_distances = mean(jump_distances);
%             q_estimate = [mean_fp_time; mean_num_jumps; mean_jump_distances];
%         end
%     end
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
    time_vec = (0.05:0.05:0.15)*0.5;
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

function dist = distance_metric(q1,q2,L, option_a_distance)
if option_a_distance %euclidian
    %NB currently no scaling between different statistics
    scaling = [400;500;0.5];
    dist = sum(((q1-q2)./scaling).^2);
    
else %kl divergence
    %equal weightings to each of the time points currently
    delx = 1;
    dist = 0;
    for j=1:length(q1(1,:))
        dist = dist + kldiv((0:delx:L)',q1(:,j)+eps,q2(:,j)+eps);
    end
end
end
