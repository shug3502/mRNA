function sample = ABC_APMC_testing

% theta ~U(-10,10)
% x| theta ~ 0.5 N(theta,1) + 0.5 N(theta,1/100) 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fake data not needed because we know the target distribution

%try making fake data 
theta_real = 0;
N=10000;
aux = rand(N,1);
bins = -10:0.1:10;
fake_data = (aux>0.5).*(theta_real+randn(N,1))+(aux<=0.5).*(theta_real+randn(N,1)*sqrt(1/100));
y = histc(fake_data,bins);
y = y./sum(y)/0.1;

num_sims = 2000;
alpha = 0.1;
M=500;
num_generations = 2;

%simulate from prior
theta = 10*(2*repmat(rand(num_sims,1),1,M)-1);
aux = repmat(rand(num_sims,1),1,M);
simulated_data = (aux>0.5).*(theta+randn(num_sims,M))+(aux<=0.5).*(theta+randn(num_sims,M)*sqrt(1/100));
x = histc(simulated_data',bins);
x = x./repmat(sum(x),length(bins),1)/0.1;

%calculate distance and decide which to keep
d = dist_to_target(x,repmat(y,1,num_sims));
weights = 1/num_sims*ones(num_sims,1);
to_keep = (d<=quantile(d,alpha));
theta = simulated_data(to_keep,:);
weights = weights(to_keep);

wmean = sum(theta.*repmat(weights,1,M))/sum(weights); %weighted mean
wvariance = 1/(length(weights)-1).*sum(repmat(weights,1,M).*(theta - repmat(wmean,length(weights),1)).^2); %weighted variance
sigma = 2*wvariance;  %weighted set of theta values

theta(2,:)

for k=1:num_generations
    k
    
    %store previous theta
    weights_store = weights;
    theta_store = theta;
    dist_store = d;
    
    N_current = length(weights);
    theta = zeros(num_sims,M);
    theta(1:N_current,1:M) = theta_store;
    previous_params = zeros(num_sims,M);
    previous_params(1:N_current,:) = theta_store;
    abc_weights = zeros(num_sims,1);
    abc_weights(1:N_current) = weights_store(1:N_current);
    d = zeros(num_sims,1);
    d(1:N_current) = dist_store(1:N_current);
    
    rand_store = zeros(num_sims-N_current,1);
    for i=1:(num_sims-N_current)
        %sample params from previous iteration
        u = rand(1);
        my_index = 1;
        %draw from discrete distribution with weights abc_weights
        while cumsum(weights_store(1:my_index))/sum(weights_store)<u
            my_index = my_index+1;
        end
        %previous_params(i,:) = par_params(my_index,p_indices); %store previous parameters in correct order to help measure displacement
        rand_store(i) = my_index;
    end
        %peturb previous parameters
        %theta(my_index,:) = theta_store(my_index,:) + sigma*randn(1);
        
        theta((N_current+1):num_sims,:) = theta_store(rand_store,:) + repmat(sigma,(num_sims-N_current),1).*repmat(randn(num_sims-N_current,1),1,M);
        
        %simulate data using these model parameters
        aux = repmat(rand(num_sims,1),1,M);
        simulated_data = (aux>0.5).*(theta+randn(num_sims,M))+(aux<=0.5).*(theta+randn(num_sims,M)*sqrt(1/100));
        x = histc(simulated_data',bins);
        x = x./repmat(sum(x),length(bins),1)/0.1;
        d = dist_to_target(x,repmat(y,1,num_sims));        
        
        %end    %repeat until N acceptances have been made
    %end
    for i=(N_current+1):N
        %Uniform prior
        prior = (theta(i,1)>-10).*(theta(i,1)<10)./20;
        kernel=1;
        kernel=kernel.*exp(-((theta(i,1)-theta_store(1:N_current,1)).^2)./(2*sigma(1).^2))./sigma(1);
        abc_weights(i) = prior./(sum(weights_store./sum(weights_store(1:N_current)).*kernel)/(sqrt(2*pi)));        
        if isnan(abc_weights(i)) || isinf(abc_weights(i))
            fprintf('Some weights are Nan or Inf\n');
            sum(weights_store)
            (theta(i,1)-previous_params(1:N_current,1)).^2
            -((theta(i,1)-previous_params(1:N_current,1)).^2)/(2*sigma(1)^2)
            -((theta(i,2)-previous_params(1:N_current,2)).^2)/(2*sigma(2)^2)
            -((theta(i,3)-previous_params(1:N_current,3)).^2)/(2*sigma(3)^2)
            abc_weights(i) = 0;
            %error('weights are NAN due to division by 0 in calculation of weights. Oops.');
        end
    end
    to_keep = (d <= quantile(d,alpha));
    theta = theta(to_keep,:);
    abc_weights = abc_weights(to_keep)./sum(abc_weights(to_keep));  %choose not to renormalise 
    d = d(to_keep);
    wmean = sum(theta.*repmat(abc_weights,1,M))/sum(abc_weights); %weighted mean
    %wvariance = (sum(abc_weights)/(sum(abc_weights)^2-sum(abc_weights.^2))).*sum(repmat(abc_weights,1,length(p_indices)).*(abc_theta - repmat(wmean,length(abc_weights),1)).^2); %weighted variance
    wvariance = 1/(length(abc_weights)-1).*sum(repmat(abc_weights,1,M).*(theta - repmat(wmean,length(abc_weights),1)).^2); %weighted variance
    if isnan(wvariance)
        error('weights are nan second time');
    end
    sigma = 2*wvariance;  %var(abc_theta.*repmat(abc_weights,1,3)); %weighted set of theta values
    %sigma_alt = 2*var(abc_theta);
    
end

sample = theta;
 
%plot the results
figure;
 subplot(2,1,1)
 hist(simulated_data,20);
subplot(2,1,2)
Num_in_bins=histc(sample,bins);
bar(bins,Num_in_bins/sum(Num_in_bins)/0.1,'histc');
hold on
vec = linspace(-10,10,201);
z = 0.5*normpdf(vec,theta_real,1) + 0.5*normpdf(vec,theta_real,1/sqrt(100));
plot(vec,z,'g--','linewidth',2);

% figure;
% bar(bins,y,'histc');
% hold on
% plot(vec,z,'g--','linewidth',2);
end

function d = dist_to_target(x,y)
%theta_real = 0;
% vec = linspace(-10,10,201);
% y = 0.5*normpdf(vec,theta_real,1) + 0.5*normpdf(vec,theta_real,1/sqrt(100));
d = sum((x-y).^2);  %(x-theta_real).^2;
end
