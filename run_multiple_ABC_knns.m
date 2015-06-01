function run_multiple_ABC_knns(num_runs)
%note that prior_params are not exactly as in ABC_knn code, only varied
%parameters
prior_params = [1.16, 0.8, 0.11, 0.42, 0.84, 0.58, 0.5, 0];
real_params = prior_params; %if real params used are different then MUST change this!
range = 0.8;
indices = [1,4,6];
posterior = zeros(11);

%num_runs = 20;
for b=1:num_runs
    [abc_theta,abc_weights]=ABC_APMC(b);
    posterior_temp = multi_dim_bin_posterior(abc_theta,abc_weights,prior_params,real_params,range,indices,0); %edited from bin_posterior
    posterior = posterior+posterior_temp;
end
posterior = posterior/num_runs;
figure;
M=11;
x1 = linspace(prior_params(indices(1))-range/2,prior_params(indices(1))+range/2,M);
x2 = linspace(prior_params(indices(2))-range/2,prior_params(indices(2))+range/2,M);
pcolor(x1,x2,posterior);
hold on 
plot(prior_params(indices(1)),prior_params(indices(2)),'ko','MarkerSize',12);
xlabel('\nu_1');
ylabel('\phi');
figure;
contourf(x1,x2,posterior);
hold on 
plot(prior_params(indices(1)),prior_params(indices(2)),'ko','MarkerSize',12);
xlabel('\nu_1');
ylabel('\phi');
end