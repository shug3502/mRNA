function mixture_model_estimation

%% fit a mixture model
%% see kamary et al 2014

%%first test hypothesis about how large the ring canal is: 10um vs 1um wide
%%fit aplha f1(x|theta) + (1-alpha) f2(x|theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=10; %repeats of the data
dt = 0.5; %1.5*10^-4;
t_end = 5; %4.5*10^-3;
%set up models to compare
%addpath ../
%addpath ../../Omero_data/scripts/
addpath ../../mRNA
addpath ../../mRNA/Model_comparison
%parameters for the mRNA model
    p1.with_anchoring = 1;
    p1.num_modes = 2;
    p1.nu1 = 1; %speed of RNP complex under active transport [zimyanin et al 2008]
    p1.nu2 = 0.5; %ratio between speed for active transport vs diffusion [zimyanin et al 2008]
    p1.lambda_1=0;   %1/0.13; %transition rate =7.69 [zimyanin et al 2008]
    p1.lambda_2 = 0.1;
    p1.omega_1= 0.4;    %1/6*(num_modes>1); %rate of falling off the microtubule [zimyanin et al 2008] $
    p1.omega_2 = 0.4;
    p1.phi = 0.58; %percentage of microtubules in posterior direction for biased angle distn [parton et $
    p1.Lx = 52; %length of cell in x direction
    p1.Ly = 37; %in y direction
    p1.Lz = 37; %in z direction
    p1.nuc_radius = 10; %radius of nucleus
    p1.theta_0 = 0; %initial angle is 0
    p1.rc_width = 1;
    p1.ellipsoid_boundary = 1;
    p1.sample_times = 0:dt:t_end; %roughly the same as time series for thresholded image
    p1.num_hist_bins=21;
    p1.zslice = [-1,1];
    p1.path_summary_stat=0;
    p1.repeats=1;
    
    p2=p1;
    %p2.ellipsoid_boundary = 0;
    %p2.phi=0.5;
    p2.rc_width=10; %compare to a much bigger ring canal
    
    %coords = get_coords('thresholded');
    xi = randn(3,10); %78 ish particles in first frame of thresholded image
    lmb = sqrt(sum(xi.^2,1));
    coords = p1.nuc_radius*xi./repmat(lmb,3,1); %generate random points at uniform on the unit sphere    
    repeats = p1.repeats; %repeats of the data
    
th_fixed = [1,0.5,0.58,0.4,0.4,0.1];    
params.model1 = @(theta) mRNA_repeats_wrapper(coords,th_fixed,p1);  %theta(1) + theta(2)*rand(k,1);
params.model2 = @(theta) mRNA_repeats_wrapper(coords,th_fixed,p2);  %theta(1) + theta(2)*randn(k,1);
params.num_sum_stats = 231; %651;
params.sigma = 0.1;
params.threshold = 0.1;
params.N=5*10^3; %number of samples
params.theta_real=[];
params.save_name = 'model_comparison_mrna_RCwidth_v2_tend5_N5000';

%data to compare to (synthetic)
tic; y = params.model1(params.theta_real); toc 
figure; plot(y(1,:,1));
tic; y2 = params.model2(params.theta_real); toc 
figure; plot(y2(1,:,1));
if params.num_sum_stats~=size(y,2)
    error('wrong number of summary stats');
end
%y = reshape(y,1,params.num_sum_stats,[]);
num_params = numel(params.theta_real);

%initialise
alpha = zeros(1,params.N); %mixing parameter
theta = zeros(num_params,params.N); %parameters of models
[alpha(1),theta(:,1)] = my_prior(num_params,1);
x = zeros(params.N,params.num_sum_stats,repeats);

parfor j=2:params.N
    if mod(j,100)==0
        fprintf('Done so far %d iterations .... \n',j);
    end
    %perturb parameters
    %alpha(j) = my_kernel(alpha(j-1),params);
    [alpha(j),theta(:,j)] = my_prior(num_params,1);
    %simulate data
    data1 = params.model1(theta(:,j));
    data2 = params.model2(theta(:,j));
    x(j,:,:) = reshape(alpha(j)*data1 + (1-alpha(j))*data2,1,params.num_sum_stats,[]);
end
%compare distances
dist_store = my_dist(repmat(y,params.N,1,1),x,ones(1,size(x,2)));

%perhaps also discard burn in?

M = ceil(params.threshold*params.N);
[~,sort_ind] = sort(dist_store);
selected = sort_ind(1:M);
alpha_abc_sample = alpha(selected);
theta_abc_sample = theta(:,selected);

%plot results and save
figure;
subplot(1+num_params,1,1);
set(gca,'fontsize',24);
plot(alpha_abc_sample,'k','linewidth',3);
ylabel('\alpha');
xlabel('t');
for i=1:num_params
    subplot(1+num_params,1,1+i);
    plot(theta_abc_sample(i,:),'k','linewidth',3);
    ylabel('\theta');
    xlabel('t');
end
print(sprintf('%s_posterior',params.save_name),'-depsc');
figure;
subplot(1+num_params,1,1);
hist(alpha_abc_sample);
xlabel('\alpha');
ylabel('frequency');
for i=1:num_params
    subplot(1+num_params,1,1+i);
    hist(theta_abc_sample(i,:));
    xlabel('\theta');
    ylabel('frequency');
end
print(sprintf('%s_trace',params.save_name),'-depsc');

save(params.save_name,'p1','p2','alpha','dist_store');

function q = my_kernel(alph,params)
%perturbation kernel
%need still inside [0,1]
q = min(max(0,alph + params.sigma*randn(1)),1);

function d = my_dist(y1,y2,weights)
%now assume y1 is Nxnum_paramsxrepeats
%weighted euclidean distance. Could probably use dist or an inbuilt fn
d = sum(sum((repmat(weights,size(y1,1),1,size(y1,3)).*bsxfun(@minus,y1,y2)).^2,2),3);

function [th1,th2] = my_prior(num_params,n)
%sample from prior
a = 0.05; %this tries to force decision between models. alpha is near 0 or 1
th1 = betarnd(a,a,1,n);
th2 = 10.^(2*rand(num_params,n)-1);
