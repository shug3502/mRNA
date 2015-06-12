function [abc_theta,abc_weights,entropy,quality,contained_in_pred_interval] = ABC_Weights_depend_on_dist(N,option_a,my_seed,data_to_read)

%created 1/6/15
%last edit 1/6/15
%AS simple as ABC
%Based on ABC_APMC
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

[params, real_params, accepted_proportion, p_accept_min, ...
    q_estimate_fake, prior_params, p_indices, prior_sigma, par_params] = initialise_parameters(data_to_read,option_a);

abc_theta = zeros(N,length(p_indices));
abc_weights = zeros(N,1);
abc_dist = zeros(N,1);
fprintf('Generation 1 begins\n');
num_generations = 1;
for i=1:N
    %initialise greater than tolerance
    %Uniform prior
    rr = rand(1,length(p_indices));
    
    %simulate parameters from the prior
    par_params(p_indices) = prior_params(p_indices)+prior_sigma(p_indices).*(rr-0.5);
    
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
% ma = max(abc_dist);
% mi = min(abc_dist);
wmean = sum(abc_theta.*repmat(abc_weights,1,length(p_indices)))/sum(abc_weights); %weighted mean
%wvariance = (sum(abc_weights)/(sum(abc_weights)^2-sum(abc_weights.^2))).*sum(repmat(abc_weights,1,length(p_indices)).*(abc_theta - repmat(wmean,length(abc_weights),1)).^2); %weighted variance
wvariance = 1/(length(abc_weights)-1).*sum(repmat(abc_weights,1,length(p_indices)).*(abc_theta - repmat(wmean,length(abc_weights),1)).^2); %weighted variance
sigma = 2*wvariance;  %var(abc_theta.*repmat(abc_weights,1,3)); %weighted set of theta values
%sigma_alt = 2*var(abc_theta);
p_accept = 1; %initialise

% figure(my_seed+1);
% subplot(3,1,1);
% plot(abc_theta(:,1),abc_theta(:,2),'o');
% hold all
% grid on
% plot(real_params(p_indices(1)),real_params(p_indices(2)),'rx','MarkerSize',12);
% set(gca, 'fontsize',14);
% xlabel('param1');
% ylabel('param2');
% subplot(3,1,2);
% plot(abc_theta(:,1),abc_theta(:,3),'o');
% hold all
% grid on
% plot(real_params(p_indices(1)),real_params(p_indices(3)),'rx','MarkerSize',12);
% set(gca, 'fontsize',14);
% xlabel('param1');
% ylabel('param3');
% subplot(3,1,3);
% plot(abc_theta(:,2),abc_theta(:,3),'o');
% hold all
% grid on
% plot(real_params(p_indices(2)),real_params(p_indices(3)),'rx','MarkerSize',12);
% set(gca, 'fontsize',14);
% xlabel('param2');
% ylabel('param3');


%First generation for t=1 is done
%Now loop over generations
while p_accept >p_accept_min
    num_generations = num_generations+1;
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
        q_estimate_candidate = summary_statistic_calculator(par_params(my_index,:),20,1,option_a);
        abc_dist(i) = distance_metric(q_estimate_candidate,q_estimate_fake,params.Lx,option_a);
        
        
        %end    %repeat until N acceptances have been made
        abc_theta(i,:) = par_params(my_index,p_indices);
    end
    for i=(N_current+1):N
        %Uniform prior
        prior = 1;
        for jj=1:length(p_indices)
            prior = prior.*(abc_theta(i,jj)>prior_params(p_indices(jj))-0.5*prior_sigma(jj)).*(abc_theta(i,jj)<prior_params(p_indices(jj))+0.5*prior_sigma(jj))./(prior_sigma(jj));
        end
        abc_weights(i) = prior./abc_dist(i).*weights_store(my_index)/sum(weights_store);
        if isnan(abc_weights(i)) || isinf(abc_weights(i))
            fprintf('Some weights are Nan or Inf\n');
            abc_weights(i) = 0;
            %error('weights are NAN due to division by 0 in calculation of weights. Oops.');
        end
    end
    to_keep = (abc_dist <= quantile(abc_dist,accepted_proportion));
    abc_theta = abc_theta(to_keep,:);
    abc_weights = abc_weights(to_keep)./sum(abc_weights(to_keep));  %choose not to renormalise
    abc_dist = abc_dist(to_keep);
%     ma = max(abc_dist);
%     mi = min(abc_dist);
    wmean = sum(abc_theta.*repmat(abc_weights,1,length(p_indices)))/sum(abc_weights); %weighted mean
    %wvariance = (sum(abc_weights)/(sum(abc_weights)^2-sum(abc_weights.^2))).*sum(repmat(abc_weights,1,length(p_indices)).*(abc_theta - repmat(wmean,length(abc_weights),1)).^2); %weighted variance
    wvariance = 1/(length(abc_weights)-1).*sum(repmat(abc_weights,1,length(p_indices)).*(abc_theta - repmat(wmean,length(abc_weights),1)).^2); %weighted variance
    if isnan(wvariance)
        error('weights are nan second time');
    end
    sigma = 2*wvariance;  %var(abc_theta.*repmat(abc_weights,1,3)); %weighted set of theta values
    %sigma_alt = 2*var(abc_theta);
    p_accept = 1/(N-N_current)*sum(to_keep((N_current+1):N));
    entropy = calculate_entropy(abc_theta,prior_params,prior_sigma,p_indices);
    
%     figure(my_seed+1);
%     subplot(3,1,1);
%     plot(abc_theta(:,1),abc_theta(:,2),'o');
%     hold all
%     grid on
%     plot(real_params(p_indices(1)),real_params(p_indices(2)),'rx','MarkerSize',12);
%     set(gca, 'fontsize',14);
%     xlabel('param1');
%     ylabel('param2');
%     subplot(3,1,2);
%     plot(abc_theta(:,1),abc_theta(:,3),'o');
%     hold all
%     grid on
%     plot(real_params(p_indices(1)),real_params(p_indices(3)),'rx','MarkerSize',12);
%     set(gca, 'fontsize',14);
%     xlabel('param1');
%     ylabel('param3');
%     subplot(3,1,3);
%     plot(abc_theta(:,2),abc_theta(:,3),'o');
%     hold all
%     grid on
%     plot(real_params(p_indices(2)),real_params(p_indices(3)),'rx','MarkerSize',12);
%     set(gca, 'fontsize',14);
%     xlabel('param2');
%     ylabel('param3');
    
end

%get rid of values outside of support of prior
abc_theta = abc_theta((abc_weights>0),:); %if weight is 0 then get rid of that parameter
%entropy gives measure of difference from uniform distn
entropy = calculate_entropy(abc_theta,prior_params,prior_sigma,p_indices);


% figure(my_seed+2);
% subplot(3,1,1);
% plot(abc_theta(:,1),abc_theta(:,2),'o');
% hold all
% grid on
% plot(real_params(p_indices(1)),real_params(p_indices(2)),'rx','MarkerSize',12);
% set(gca, 'fontsize',14);
% xlabel('param1');
% ylabel('param2');
% subplot(3,1,2);
% plot(abc_theta(:,1),abc_theta(:,3),'o');
% hold all
% grid on
% plot(real_params(p_indices(1)),real_params(p_indices(3)),'rx','MarkerSize',12);
% set(gca, 'fontsize',14);
% xlabel('param1');
% ylabel('param3');
% subplot(3,1,3);
% plot(abc_theta(:,2),abc_theta(:,3),'o');
% hold all
% grid on
% plot(real_params(p_indices(2)),real_params(p_indices(3)),'rx','MarkerSize',12);
% set(gca, 'fontsize',14);
% xlabel('param2');
% ylabel('param3');

%post-processing to check quality of posterior
[entropy, quality, contained_in_pred_interval]= post_processing(abc_theta,abc_weights,prior_params,real_params,prior_sigma,p_indices);

fname = sprintf('ABC_APMC_output%d.txt',my_seed);
fileID = fopen(fname,'w');
fprintf(fileID,'%f \n',abc_theta);
fclose('all');

toc
end
