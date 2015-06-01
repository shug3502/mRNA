function posterior = multi_dim_bin_posterior(abc_theta,prior_params,range,indices)
%prior params should be a vector with the centres of the uniform
%distribution for the prior
%range is the spread of the normal prior. Assumed this is the same for all
%the parameters
%computes multidimensional histogram

if nargin~=4
    fprintf('Using defaults \n');
    fname = sprintf('ABC_APMC_output2.txt');   %('Simplest_ABC_with_moments_output2.txt');
    fileID = fopen(fname,'r');
    abc_theta = fscanf(fileID, '%f');
    N = length(abc_theta)/3;
    abc_theta = reshape(abc_theta,N,3);
    fclose('all');
    prior_params = [1.16, 0.8, 0.11, 0.42, 0.84, 0.58, 0.5, 0];
    range = 0.4;
    indices = [1,4,6];
end
N = length(abc_theta(:,1));
M=11;
abc_theta
%split into bins
edges_vecs = [linspace(prior_params(indices(1))-range/2,prior_params(indices(1))+range/2,M);
    linspace(prior_params(indices(2))-range/2,prior_params(indices(2))+range/2,M);
    linspace(prior_params(indices(3))-range/2,prior_params(indices(3))+range/2,M)];

[count, edges, mid, loc] = histcn(abc_theta, edges_vecs(1,:),edges_vecs(2,:),edges_vecs(3,:));
posterior = count/N; 

figure;
imagesc(mid{1:2},posterior(:,:,ceil(end/2)));

figure;
pcolor(0.5*(edges_vecs(1,1:M-1)+edges_vecs(1,2:M)),0.5*(edges_vecs(2,1:M-1)+edges_vecs(2,1:M-1)), mean(posterior,3));

figure;
contourf(0.5*(edges_vecs(1,1:M-1)+edges_vecs(1,2:M)),0.5*(edges_vecs(2,1:M-1)+edges_vecs(2,1:M-1)), mean(posterior,3));

% figure;
% contourf(posterior)
% figure;
% pcolor(posterior)

end