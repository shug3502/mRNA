function run_multiple_ABC_knns(num_runs)
%note that prior_params are not exactly as in ABC_knn code, only varied
%parameters
prior_params = [1.16, 0.42, 0.58];
range = 0.4;
indices = [1,3];
posterior = zeros(11);

%num_runs = 20;
parfor b=1:num_runs
    [abc_theta,~]=ABC_APMC(b);
    posterior_temp = bin_posterior(abc_theta,prior_params,range,indices);
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