function gradient_est = estimate_gradient(theta,data_summary,opts,current_cost);

%created 2/2/2016 JH
%last edit 2/2/2016 

%theta = parameters
%data_summary = summary statistics of data (synthetic or real)
%opts = options structure
%current_cost = cost at current parameters theta
%using a foward euler estimate reduces new estimates of f(X(theta)), although less accurate
%%%%%%%%%%%%%%%%%%%%%%%%%%

epsilon = (10^-4)/2*ones(1,3);  %step size divided by 2;

gradient_est = zeros(1,numel(theta));
for j=1:numel(theta)
	vec = zeros(1,numel(theta));
	vec(j) = 1;
	gradient_est(j) = (evaluate_cost_fn(theta + epsilon.*vec, data_summary, opts) - current_cost)/(epsilon*vec');
% evaluate_cost_fn(theta - epsilon.*vec, data_summary, opts))/(epsilon*vec');
end

end
