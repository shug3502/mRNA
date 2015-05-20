function [abc_theta,abc_weights] = Simplest_ABC_with_moments(my_seed)

%created 13/5/15
%last edit 20/5/15
%AS simple as ABC
%implements a version of ABC via rejection sampling to fit parameters to
%some 'data' which is first simulated from the model
%Model used is velocity jump process with 1 mode
%see Turner 2012 for ABC population monte carlo tutorial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic; 
rng(my_seed);
close all

%Fake parameters
params.nu1 = 1.0;
params.nu2 = 0;
params.lambda = 7.7;
params.omega = 0;
params.phi = 0.58; 
params.theta_0 = 0; 
params.x_0=0.5;

L=30;
delx=1;
real_params = [params.nu1, params.nu2, params.lambda, params.omega, params.phi, params.theta_0, params.x_0];

%Generate fake data
%calculate appropriate summary statistic - we choose MFPT
q_estimate_fake = summary_statistic_calculator(params,100,0,65)

%Choose tolerance sequence
num_generations = 3;
%delta values chosen by looking at the standard deviations of the summary
%statistics as parameters vary
delta = [10,5,2,0.5,0.1]; %[40,20,10,5];

%At t=1 for first generation
N=100; 
%N is # particles at each generation
%Note that algorithm does not perform well when N is small. Needs to pick
%sensible parameters, or acceptance rate is very very low. 

%create while loop

%set prior
prior_params = [0.9, 0, 7.7, 0, 0.58, 0, 0.5];
prior_sigma = [0.4, 0.4, 0.4]; %sd of gaussian or spread around mean of uniform
p_indices = [1, 3, 5];

n=0;
abc_theta = zeros(N,length(p_indices));
abc_weights = zeros(N,1);
fprintf('Generation 1 begins\n');
for i=1:N
    %initialise greater than tolerance
    q_estimate_candidate =[ones(1,3); zeros(L/delx,3)]; % ones(1+L/delx,3)/(1+L/delx); %
    %kldiv((0:delx:L)',q_estimate_candidate(:,1)+eps,q_estimate_fake(:,1)+eps)
    par_params = prior_params;
while kldiv((0:delx:L)',q_estimate_candidate(:,1)+eps,q_estimate_fake(:,1)+eps) > delta(1)    %compare the summary statistic to that from original data

    n=n+1;
%Gaussian prior
%rr = randn(1,3);
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
    q_estimate_candidate = summary_statistic_calculator(par_params,20,1,4);
    if n>10^5
        fprintf('Thats more than enough parameters to look at for now.\n')
        break
    end
end    %repeat until N acceptances have been made
abc_theta(i,:) = par_params(p_indices);
abc_weights(i) = 1/N;
end
sigma = 2*var(abc_theta);


figure(1);
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
for tau=2:num_generations
    %store previous theta
    theta_store = abc_theta;
    weights_store = abc_weights;
    par_params = repmat(prior_params,N,1);
    par_params(:,p_indices) = theta_store;
    n=0;
    tau
for i=1:N
    %initialise greater than tolerance
    q_estimate_candidate =[ones(1,3); zeros(L/delx,3)];
while kldiv((0:delx:L)',q_estimate_candidate(:,1)+eps,q_estimate_fake(:,1)+eps)>delta(tau)    %compare the summary statistic to that from original data
    n=n+1

    %simulate parameters from the prior    
%     candidate.nu1 = 0.4 + 0.1*rr(1,j);
%     candidate.nu2 = 0.08 + 0.02*rr(2,j);
%     candidate.lambda = 7.7; % + 2*randn(1);
%     candidate.omega = 1/6; %and this  %+ 0.5*randn(1);
%     candidate.phi = 0.58; %fix this too %+0.2*randn(1);
%     candidate.theta_0 = 0; %these are fixed
%     candidate.x_0 = 0.5; %these are fixed

%sample params from previous iteration
u = rand(1);
my_index = 1;
%draw from discrete distribution with weights abc_weights
while cumsum(weights_store(1:my_index))/sum(weights_store)<u
    my_index = my_index+1;
end


%peturb previous parameters
par_params(my_index,p_indices) = par_params(my_index,p_indices) + sigma.*randn(1,length(p_indices)); 

    %simulate data using these model parameters
    %Calculate summary statistic (MFPT)
    q_estimate_candidate = summary_statistic_calculator(par_params(my_index,:),20,1,4);
    if n>10^5
        fprintf('Thats more than enough parameters to look at for now.\n')
        break
    end
 if abs(sum(q_estimate_candidate(:,1)+eps)-1)>10^-4
q_estimate_candidate     
toc;
 end
end    %repeat until N acceptances have been made
abc_theta(i,:) = par_params(my_index,p_indices);
%Gaussian prior
% prior = exp(-((abc_theta(i,1)-prior_params(p_indices(1))^2)/(2*prior_sigma(1)^2)))*...
%     exp(-((abc_theta(i,2)-prior_params(p_indices(2)))^2)/(2*prior_sigma(2)^2))*...
%     exp(-((abc_theta(i,3)-prior_params(p_indices(3)))^2)/(2*prior_sigma(3)^2));
%Uniform prior
prior = 1;
for jj=1:3
prior = prior.*(abc_theta(i,jj)>prior_params(p_indices(jj))-0.5*prior_sigma(jj)).*(abc_theta(i,jj)<prior_params(p_indices(jj))+0.5*prior_sigma(jj))./(prior_sigma(jj));
end
%fixed dependence on current weights to depend on previous weights
abc_weights(i) = prior./(sum(weights_store.*exp(-((abc_theta(:,1)-theta_store(:,1)).^2)/(2*sigma(1)^2))...
    .*exp(-((abc_theta(:,2)-theta_store(:,2)).^2)/(2*sigma(2)^2))...
    .*exp(-((abc_theta(:,3)-theta_store(:,3)).^2)/(2*sigma(3)^2)))/(sigma(1)*sigma(2)*sigma(3))^2);
if isnan(abc_weights)
    abc_weights
    error('weights are NAN due to division by 0 in calculation of weights. Oops.');
end
end
  sigma = 2*var(abc_theta);  
    

figure(1);
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
length(abc_theta)
%get rid of values outside of cupport of prior
abc_theta = abc_theta((abc_weights>0),:); %if weight is 0 then get rid of that parameter
length(abc_theta)


figure(2);
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

fname = sprintf('Simplest_ABC_with_moments_output%d.txt',my_seed);
fileID = fopen(fname,'w');
fprintf(fileID,'%f \n',abc_theta);
fclose('all');

toc
end

function [q_estimate] = summary_statistic_calculator(par_params,N,is_parallel, r_random)
%runs velocityjump2D_ModesInput which is a velocity jump process for a single particle
%runs this many times for many particles and returns simulated estimates of
%mfpt

delx = 1; %bin size;
if is_parallel
    params.nu1 = par_params(1);
    params.nu2 = par_params(2);
    params.lambda = par_params(3);
    params.omega = par_params(4);
    params.phi = par_params(5); 
    params.theta_0 = par_params(6);     
    params.x_0 = par_params(7); 
else
    params = par_params;
end

L=30;
time_vec = (0.05:0.05:0.05);
t=time_vec*60^2;
t_max = 1;  
l_t = length(time_vec);
jumps = zeros(N,10^4);
q_estimate = zeros(L/delx+1,l_t);
xpos = zeros(N,10^4);
xpos_discrete_time = zeros(N,10^4);

for j=1:N
    [~, ~, ~, ~, xpos_temp, jump_temp] = velocityjump2D_ModesInput(max(time_vec), params, 1, 0, 1, 0, r_random*j);
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
        ename = sprintf('Simplest_ABC_with_moments_errors.txt');
fileIDerror = fopen(ename,'w');
fprintf(fileIDerror,'%f \n',xpos_discrete_time);
fclose('all');
        error('outside of 0:L');
    end
    q_estimate(:,w) = Num_in_bins/N/delx; %estimate of q at time T
%     mean_position(w) = (delx/2:delx:(L-delx/2))*q_estimate(1:(end-1),w)*delx+delx*L*q_estimate(end,w); %use centre of each bin
%     second_moment(w) = (delx/2:delx:(L-delx/2)).^2*q_estimate(1:(end-1),w)*delx+delx*L^2*q_estimate(end,w); %use centre of each bin
end

end