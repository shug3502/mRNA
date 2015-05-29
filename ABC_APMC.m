function [abc_theta,abc_weights] = ABC_APMC(my_seed)

%created 20/5/15
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
rng(my_seed);
close all

%Fake parameters
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

real_params = [params.nu1, params.nu2, params.lambda_1, params.lambda_2, params.omega_1, params.omega_2, params.phi, params.x_0, params.Lx, params.Ly, params.nuc_radius, params.theta_0];

%Generate fake data
%calculate appropriate summary statistic - we choose MFPT
q_estimate_fake = summary_statistic_calculator(params,1000,0)

%Choose tolerance sequence
accepted_proportion = 0.5; %alpha
%At t=1 for first generation
N=2000;

p_accept_min = 0.1; % 1%

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
num_generations = 1;
for i=1:N
    %initialise greater than tolerance
    %Uniform prior
    rr = rand(1,3);
    
    %simulate parameters from the prior
    %     candidate.nu1 = 0.4 + 0.1*rr(1,j);
    %     candidate.nu2 = 0.08 + 0.02*rr(2,j);
    %     candidate.lambda = 7.7; % + 2*randn(1);
    %     candidate.omega = 1/6; %and this  %+ 0.5*randn(1);
    %     candidate.phi = 0.58; %fix this too %+0.2*randn(1);
    %     candidate.theta_0 = 0; %these are fixed
    %     candidate.x_0 = 0.5; %these are fixed
    par_params(p_indices) = prior_params(p_indices)+prior_sigma.*(rr-0.5);
    
    %simulate data using these model parameters
    %Calculate summary statistic (MFPT)
    q_estimate_candidate = summary_statistic_calculator(par_params,20,1);
    abc_dist(i) = distance_metric(q_estimate_candidate,q_estimate_fake,params.Lx); %distance of proposed S(x) from S(x_obs)
    %end    %repeat until N acceptances have been made
    
    abc_theta(i,:) = par_params(p_indices);
    abc_weights(i) = 1/N;
end

%keep now parameters with distances less than the accepted quantile
to_keep = (abc_dist <= quantile(abc_dist,accepted_proportion));
abc_theta = abc_theta(to_keep,:);
abc_weights = abc_weights(to_keep)./sum(abc_weights(to_keep));
abc_dist = abc_dist(to_keep);
ma = max(abc_dist);
mi = min(abc_dist);
wmean = sum(abc_theta.*repmat(abc_weights,1,length(p_indices)))/sum(abc_weights); %weighted mean
%wvariance = (sum(abc_weights)/(sum(abc_weights)^2-sum(abc_weights.^2))).*sum(repmat(abc_weights,1,length(p_indices)).*(abc_theta - repmat(wmean,length(abc_weights),1)).^2); %weighted variance
wvariance = 1/(length(abc_weights)-1).*sum(repmat(abc_weights,1,length(p_indices)).*(abc_theta - repmat(wmean,length(abc_weights),1)).^2); %weighted variance
sigma = 2*wvariance;  %var(abc_theta.*repmat(abc_weights,1,3)); %weighted set of theta values
%sigma_alt = 2*var(abc_theta);
p_accept = 1; %initialise

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


%First generation for t=1 is done
%Now loop over generations
while p_accept >p_accept_min
    num_generations = num_generations+1
    %store previous theta
    weights_store = abc_weights;
    theta_store = abc_theta;
    dist_store = abc_dist;
    
    N_current = length(abc_weights);
    par_params = repmat(prior_params,N,1);
    par_params(1:N_current,p_indices) = theta_store;
    abc_theta = zeros(N,length(p_indices));
    abc_theta(1:N_current,1:length(p_indices)) = theta_store;
    previous_params = zeros(N,length(p_indices));
    previous_params(1:N_current,:) = theta_store;
    abc_weights = zeros(N,1);
    abc_weights(1:N_current) = weights_store(1:N_current);
    abc_dist = zeros(N,1);
    abc_dist(1:N_current) = dist_store(1:N_current);
    
    for i=(N_current+1):N
        %sample params from previous iteration
        u = rand(1);
        my_index = 1;
        %draw from discrete distribution with weights abc_weights
        while cumsum(weights_store(1:my_index))/sum(weights_store)<u
            my_index = my_index+1;
        end
        previous_params(i,:) = par_params(my_index,p_indices); %store previous parameters in correct order to help measure displacement
        
        %peturb previous parameters
        par_params(my_index,p_indices) = par_params(my_index,p_indices) + sigma.*randn(1,length(p_indices));
        
        %simulate data using these model parameters
        %Calculate summary statistic (MFPT)
        q_estimate_candidate = summary_statistic_calculator(par_params(my_index,:),20,1);
        abc_dist(i) = distance_metric(q_estimate_candidate,q_estimate_fake,params.Lx);
        
        
        %end    %repeat until N acceptances have been made
        abc_theta(i,:) = par_params(my_index,p_indices);
    end
    for i=(N_current+1):N
        %Uniform prior
        prior = 1;
        for jj=1:3
            prior = prior.*(abc_theta(i,jj)>prior_params(p_indices(jj))-0.5*prior_sigma(jj)).*(abc_theta(i,jj)<prior_params(p_indices(jj))+0.5*prior_sigma(jj))./(prior_sigma(jj));
        end
        abc_weights(i) = prior./(sum(weights_store./sum(weights_store).*exp(-((abc_theta(i,1)-previous_params(1:N_current,1)).^2)/(2*sigma(1)^2))...
            .*exp(-((abc_theta(i,2)-previous_params(1:N_current,2)).^2)/(2*sigma(2)^2))...
            .*exp(-((abc_theta(i,3)-previous_params(1:N_current,3)).^2)/(2*sigma(3)^2)))/(sqrt(2*pi)^length(p_indices)*(sigma(1)*sigma(2)*sigma(3))^2));
        if isnan(abc_weights(i)) || isinf(abc_weights(i))
            fprintf('Some weights are Nan or Inf\n');
            sum(weights_store)
            -((abc_theta(i,1)-previous_params(1:N_current,1)).^2)/(2*sigma(1)^2)
            -((abc_theta(i,2)-previous_params(1:N_current,2)).^2)/(2*sigma(2)^2)
            -((abc_theta(i,3)-previous_params(1:N_current,3)).^2)/(2*sigma(3)^2)
%             abc_weights(i) = exp(log(prior)-log((sum(weights_store./sum(weights_store).*exp(-((abc_theta(i,1)-previous_params(1:N_current,1)).^2)/(2*sigma(1)^2))...
%             .*exp(-((abc_theta(i,2)-previous_params(1:N_current,2)).^2)/(2*sigma(2)^2))...
%             .*exp(-((abc_theta(i,3)-previous_params(1:N_current,3)).^2)/(2*sigma(3)^2)))/(sqrt(2*pi)^length(p_indices)*(sigma(1)*sigma(2)*sigma(3))^2))))
            abc_weights(i) = 0;
            %error('weights are NAN due to division by 0 in calculation of weights. Oops.');
        end
    end
    to_keep = (abc_dist <= quantile(abc_dist,accepted_proportion));
    abc_theta = abc_theta(to_keep,:);
    abc_weights = abc_weights(to_keep)./sum(abc_weights(to_keep));
    abc_dist = abc_dist(to_keep);
    ma = max(abc_dist);
    mi = min(abc_dist);
    wmean = sum(abc_theta.*repmat(abc_weights,1,length(p_indices)))/sum(abc_weights); %weighted mean
    %wvariance = (sum(abc_weights)/(sum(abc_weights)^2-sum(abc_weights.^2))).*sum(repmat(abc_weights,1,length(p_indices)).*(abc_theta - repmat(wmean,length(abc_weights),1)).^2); %weighted variance
    wvariance = 1/(length(abc_weights)-1).*sum(repmat(abc_weights,1,length(p_indices)).*(abc_theta - repmat(wmean,length(abc_weights),1)).^2); %weighted variance
    if isnan(wvariance)
        error('weights are nan second time');
    end
    sigma = 2*wvariance;  %var(abc_theta.*repmat(abc_weights,1,3)); %weighted set of theta values
    %sigma_alt = 2*var(abc_theta);
    p_accept = 1/(N-N_current)*sum(to_keep((N_current+1):N)); %set this
    entropy = calculate_entropy(abc_theta,prior_params,prior_sigma,p_indices);
    
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
    
end

%length(abc_theta)
%get rid of values outside of support of prior
abc_theta = abc_theta((abc_weights>0),:); %if weight is 0 then get rid of that parameter
%length(abc_theta)
entropy = calculate_entropy(abc_theta,prior_params,prior_sigma,p_indices);


figure(my_seed+2);
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

fname = sprintf('ABC_APMC_output%d.txt',my_seed);
fileID = fopen(fname,'w');
fprintf(fileID,'%f \n',abc_theta);
fclose('all');

toc
end

function [q_estimate] = summary_statistic_calculator(par_params,N,is_parallel)
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

L=params.Lx;
time_vec = (0.05:0.05:0.15)*0.2;
t=time_vec*60^2;
l_t = length(time_vec);
jumps = zeros(N,10^4);
q_estimate = zeros(L/delx+1,l_t);
xpos = zeros(N,10^4);
xpos_discrete_time = zeros(N,10^4);

for j=1:N
    [~, ~, ~, ~, xpos_temp, ~, jump_temp] = velocityjump2D_with_nucleus(max(time_vec), params, 1, 1, 0);
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
function dist = distance_metric(q1,q2,L)
%equal weightings to each of the time points currently
delx = 1;
dist = 0;
for j=1:length(q1(1,:))
    dist = dist + kldiv((0:delx:L)',q1(:,j)+eps,q2(:,j)+eps);
end
end
