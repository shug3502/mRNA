function [entropy, quality, contained_in_pred_interval]= post_processing(abc_theta,abc_weights,prior_params,real_params,prior_sigma,p_indices)

%entropy gives measure of difference from uniform distn
entropy = calculate_entropy(abc_theta,prior_params,prior_sigma,p_indices);

%post-processing to check quality of posterior
[~, quality] = multi_dim_bin_posterior(abc_theta,abc_weights,prior_params,real_params,prior_sigma(1),p_indices,1);
%find 95% interval for posterior
box = zeros(length(p_indices),2);
contained_in_pred_interval = 1;
for k=1:3
    box(k,:) = [quantile(abc_theta(:,k),0.025),quantile(abc_theta(:,k),1-0.025)];
    contained_in_pred_interval = contained_in_pred_interval*(real_params(p_indices(k))>box(k,1))*(real_params(p_indices(k))<box(k,2));
end


end