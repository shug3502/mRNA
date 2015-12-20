function wrapper_harrison_GP()
%script to wrap GP function

%train for output m=1
m=1:7;

%extract data 
%M = csvread('train.csv'); %import data

x = 10*rand(1000,7);
z = sin(abs(x)) + sqrt(abs(x));
figure; plot(x,z,'o');
M = [z,x];
l_theta = 7; %number of output params
n = size(M,1);
train_prop = 0.8; % how much data to train on
Xtrain = M(1:train_prop*n,(l_theta+1):end);
Xtest = M((train_prop*n+1):n,(l_theta+1):end);
Y = M(1:train_prop*n,m);
Ytest = M((train_prop*n+1):n,m);

%specify covariance function
%squared exponential kernel
params.sigma = 0.1; %sigma squared
params.l = 1;
k = @(x1,x2) my_cov(x1,x2);  % params.sigma*exp(-0.5/params.l^2*(x1-x2).^2);

sn = 0.1; %variance of assumed additive noise

[fBarStar, CVfStar] = harrison_GP(Xtrain,Xtest,Y,k,sn);
size(fBarStar)

error = mean((fBarStar - Ytest).^2,1)

figure; plot(Ytest(:,1),fBarStar(:,1),'o');
end
function kk = my_cov(x1,x2)
params.sigma = 0.1; %sigma squared
params.l = 1;
aux = zeros();
for i=1:size(x1,1)
    for j=1:size(x2,1)
        aux(i,j) = sum((x1(i,:) - x2(j,:)).^2,2);
    end
end
kk =  params.sigma*exp(-0.5/params.l^2*aux);
end
