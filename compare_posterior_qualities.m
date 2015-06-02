function compare_posterior_qualities(N,num_runs)
%compare quality of posterior from ABC with ABC_APMC and ABC_weights
%dependent on distance
prior_params = [1.16, 0.8, 0.11, 0.42, 0.84, 0.58, 0.5, 0];
real_params = prior_params; %if real params used are different then MUST change this!
real_params(1) = 1;
range = 0.8;
indices = [1,4,6];

proportion_in_interval = 0;
%num_runs = 20;
for b=1:num_runs
    [~,~,entropy,quality,contained_in_pred_interval] =ABC_APMC(N,b);
    proportion_in_interval = proportion_in_interval+contained_in_pred_interval;
end
posterior = posterior/num_runs;



end