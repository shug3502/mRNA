function ABC_testing

% theta ~U(-10,10)
% x| theta ~ 0.5 N(theta,1) + 0.5 N(theta,1/100) 
% 

%fake data not needed because we know the target distribution


num_sims = 10000;

%simulate from prior
theta = 10*(2*rand(num_sims,1)-1);
aux = rand(num_sims,1);
simulated_data = (aux>0.5).*(theta+randn(num_sims,1))+(aux<=0.5).*(theta+randn(num_sims,1)*sqrt(1/100));

d = dist_to_target(x);
to_keep = ;
% figure;
% hist(simulated_data,20);


end

function d = dist_to_target(x)
theta_real = 0;
% vec = linspace(-10,10,201);
% y = 0.5*normpdf(vec,theta_real,1) + 0.5*normpdf(vec,theta_real,1/sqrt(100));
d = (x-theta_real).^2; %(x-y).^2;
end