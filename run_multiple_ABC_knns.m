function run_multiple_ABC_knns
prior_params = [1.16, 0, 0.42, 0, 0.58, 0, 0.5];
range = 0.4;
indices = [1,2];
posterior = zeros(21);
parfor b=1:4
    [abc_theta,~]=ABC_knn(b);
    posterior_temp = bin_posterior(abc_theta,prior_params,range,indices);
    posterior = posterior+posterior_temp;
end
posterior = posterior/4;
figure;
pcolor(posterior);
end