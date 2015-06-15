function logistic_growth_fitting(N)
tic;
logistic_params.r = 0.01;
logistic_params.K = 10;
[T,Y] = ode45(@logistic_ode,[0 60*3600],0.2,[],logistic_params);
figure;
plot(T,Y(:,1),'-o')


data = [0,15; 24,25; 32,35; 40,50; 46,75; 51,85; 54,115; 60,120];
data(:,1)=data(:,1)*60^2; %convert to seconds from hours
toc
tic;
rng(12);
close all

accepted_proportion=0.5;
prior_params = [0.5,100];
prior_sigma = [1,200];
real_params = prior_params;
par_params = prior_params;
p_indices=1:2;

abc_theta = zeros(N,length(p_indices));
abc_weights = zeros(N,1);
abc_dist = zeros(N,1);
fprintf('Generation 1 begins\n');

for i=1:N
    %initialise greater than tolerance
    %Uniform prior
    rr = rand(1,length(p_indices));
    
    %simulate parameters from the prior
    par_params(p_indices) = prior_params(p_indices)+prior_sigma(p_indices).*(rr-0.5)
    
    %simulate data using these model parameters
    %Calculate summary statistic (MFPT)
    q_estimate_candidate = logistic_summary_statistic_calculator(par_params,data);
    abc_dist(i) = logistic_distance_metric(q_estimate_candidate,data); %distance of proposed S(x) from S(x_obs)
    %end    %repeat until N acceptances have been made
    
    abc_theta(i,:) = par_params(p_indices);
    abc_weights(i) = 1/N;
end

%keep now parameters with distances less than the accepted quantile
to_keep = (abc_dist <= quantile(abc_dist,accepted_proportion));
abc_theta = abc_theta(to_keep,:);
abc_weights = abc_weights(to_keep)./sum(abc_weights(to_keep));
abc_dist = abc_dist(to_keep);

figure(1);
subplot(3,1,1);
plot(abc_theta(:,1),abc_theta(:,2),'o');
hold all
grid on
plot(real_params(p_indices(1)),real_params(p_indices(2)),'rx','MarkerSize',12);
set(gca, 'fontsize',14);
xlabel('param1');
ylabel('param2');

[entropy, quality, contained_in_pred_interval]= post_processing(abc_theta,abc_weights,prior_params,real_params,prior_sigma,p_indices);

toc
end

function x_summary = logistic_summary_statistic_calculator(par_params,data)
logistic_params.r = par_params(1);
logistic_params.K = par_params(2);
[T,Y] = ode45(@logistic_ode,data(:,1),data(1,2),[],logistic_params);

x_summary = Y;
end

function dist = logistic_distance_metric(x_summary, data)
dist = sqrt(sum((x_summary-data(:,1)).^2));
end

function dy = logistic_ode(t,y,logistic_params)
dy = logistic_params.r*y.*(1-y/logistic_params.K);
end