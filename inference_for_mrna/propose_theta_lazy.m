function [early_stopping, a_star, params] = propose_theta_lazy(params,opts)
%draw params
rr = rand(3,1);
theta = 10.^[2*rr(1)-1, 2*rr(2)-1, -rr(3)]; %prior uniform on log of parameters
params([1,4,6]) = theta;

av_jump = theta(1)/theta(2); %calculate average jump distance
early_stopping = (av_jump>opts.max_jump_length);
a_star = early_stopping*opts.continue_prob + (1-early_stopping);
end


