%what is the kl divergence between the prior and the posterior after lazy
%abc? Gives measure of gain in information

addpath ../
load('lazy_abcv3.mat');
num_params = 3;
for j=1:num_params
    [nbins_prior,c] = hist(theta_store(:,j));
[nbins_post] = hist(posterior(:,j),c);
A = nbins_prior/sum(nbins_prior);
B = nbins_post/sum(nbins_post);
k = kldiv(c,A+eps,B+eps);
[h, p] = kstest2(posterior(:,j),theta_store(:,j)); %test if distributions are different (although why are you doing this to test if prior and posterior are different?)
fprintf('For the %dth parameter, kl divergence is: %f and ks test returns %d at a p-value of %f \n \n', j,k,h,p);
end
