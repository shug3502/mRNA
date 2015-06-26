function sample = ABC_testing

% theta ~U(-10,10)
% x| theta ~ 0.5 N(theta,1) + 0.5 N(theta,1/100) 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fake data not needed because we know the target distribution

%try making fake data 
theta_real = 0;
N=10000;
aux = rand(N,1);
bins = -10:0.1:10;
fake_data = (aux>0.5).*(theta_real+randn(N,1))+(aux<=0.5).*(theta_real+randn(N,1)*sqrt(1/100));
y = histc(fake_data,bins);
y = y./sum(y)/0.1;

num_sims = 200000;
alpha = 0.001;
M=500;

%simulate from prior
theta = 10*(2*repmat(rand(num_sims,1),1,M)-1);
aux = repmat(rand(num_sims,1),1,M);
simulated_data = (aux>0.5).*(theta+randn(num_sims,M))+(aux<=0.5).*(theta+randn(num_sims,M)*sqrt(1/100));
x = histc(simulated_data',bins);
x = x./repmat(sum(x),length(bins),1)/0.1;

%calculate distance and decide which to keep
d = dist_to_target(x,repmat(y,1,num_sims));
to_keep = (d<=quantile(d,alpha));
sample = simulated_data(to_keep);
 
%plot the results
figure;
 subplot(2,1,1)
 hist(simulated_data,20);
subplot(2,1,2)
Num_in_bins=histc(sample,bins);
bar(bins,Num_in_bins/sum(Num_in_bins)/0.1,'histc');
hold on
vec = linspace(-10,10,201);
z = 0.5*normpdf(vec,theta_real,1) + 0.5*normpdf(vec,theta_real,1/sqrt(100));
plot(vec,z,'g--','linewidth',2);

% figure;
% bar(bins,y,'histc');
% hold on
% plot(vec,z,'g--','linewidth',2);
end

function d = dist_to_target(x,y)
%theta_real = 0;
% vec = linspace(-10,10,201);
% y = 0.5*normpdf(vec,theta_real,1) + 0.5*normpdf(vec,theta_real,1/sqrt(100));
d = sum((x-y).^2);  %(x-theta_real).^2;
end
