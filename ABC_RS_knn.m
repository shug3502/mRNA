function [abc_theta,abc_weights,entropy,quality,contained_in_pred_interval] = ABC_RS_knn(N,option_a,my_seed,data_to_read)

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

[params, real_params, accepted_proportion, p_accept_min, ...
    q_estimate_fake, prior_params, p_indices, prior_sigma, par_params] = initialise_parameters(data_to_read,option_a);


abc_theta = zeros(N,length(p_indices));
abc_weights = zeros(N,1);
abc_dist = zeros(N,1);
fprintf('Generation 1 begins\n');

for i=1:N
    if mod(i,10)==0
        i
    end
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

if length(p_indices)>=3
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

[entropy, quality, contained_in_pred_interval]= post_processing(abc_theta,abc_weights,prior_params,real_params,prior_sigma,p_indices);

% fname = sprintf('ABC_APMC_output%d.txt',my_seed);
% fileID = fopen(fname,'w');
% fprintf(fileID,'%f \n',abc_theta);
% fclose('all');

toc
end
