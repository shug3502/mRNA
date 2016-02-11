function my_loss = abc_knn(data,params)
%created 22/12/15
%last edit 15/1/16

%perform abc via knn
%%%%%%%%%%%%%%%%%%%%

%load in data in train/test format
prop = params.prop;
ntr = size(data,1)*prop; nte = size(data,1)*(1-prop);                        % number of training and test points
%my_loss = zeros(params.size_theta,1);

xtr = data(1:ntr,(1+params.size_theta):end); %train inputs
ytr = data(1:ntr,1:params.size_theta); %train responses
xte = data((1+ntr):(ntr+nte),(1+params.size_theta):end); %test inputs
yte = data((1+ntr):(ntr+nte),1:params.size_theta); %test responses

figure;
plot3(ytr(:,1),ytr(:,2),ytr(:,3),'o');
size(ytr)
check = knnsearch(ytr,[1.16,1.26,0.58],'K',10,'distance','euclidean');



mean(ytr(check,:),1)

%X and Y have observations as rows and variables as columns
idx = knnsearch(xtr,xte,'K',params.k,'distance','euclidean');
%idx now contains indices of all nearest neighbours for each observation
%these should now be used to estimate the posterior

ypred = zeros(size(yte));
for j=1:size(xte,1) %possibly better way than a loop
temp = ytr(idx(j,:),:); %matrix k by size_ss
ypred(j,:) = mean(temp,1);
end

my_loss = mean((ypred./yte - 1).^2)'

%Test for default parameters
default_params = [1.16, 0, 0, 0.42, 0, 0.58, 0];
addpath ../../
nsamples=1;
M = zeros(nsamples,1113);
for j=1:nsamples
q = reshape(summary_statistic_calculator_3D(default_params,5,1,0),1,[]);
M(j,:) = q;
end
idx2 = knnsearch(xtr,M,'K',params.k,'distance','euclidean');
for j=1:nsamples
    temp2 = ytr(idx2(j,:),:);
z(j,:) = mean(temp2,1);
end
mean(z,1)
figure;
for j=1:3
subplot(3,1,j);
hist(temp2(:,j))
end
