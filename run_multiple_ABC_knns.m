function run_multiple_ABC_knns(num_runs)
%note that prior_params are not exactly as in ABC_knn code, only varied
%parameters
prior_params = [1.16, 0.42, 0.58];
range = 0.4;
indices = [1,2];
posterior = zeros(11);

%num_runs = 20;
parfor b=1:num_runs
    [abc_theta,~]=ABC_knn(b);
    posterior_temp = bin_posterior(abc_theta,prior_params,range,indices);
    posterior = posterior+posterior_temp;
end
posterior = posterior/num_runs
figure;
pcolor(posterior);
figure;
contourf(posterior);
end