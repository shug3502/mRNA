function posterior = bin_posterior(abc_theta,prior_params,range,indices)
%prior params should be a vector with the centres of the uniform
%distribution for the prior
%range is the spread of the normal prior. Assumed this is the same for all
%the parameters

if nargin~=4
    fprintf('Using defaults');
    fname = sprintf('ABC_knn_output4.txt');   %('Simplest_ABC_with_moments_output2.txt');
    fileID = fopen(fname,'r');
    abc_theta = fscanf(fileID, '%f');
    N = length(abc_theta)/3;
    abc_theta = reshape(abc_theta,N,3);
    fclose('all');
    prior_params = [1.16,0.42,0.58];
    range = 0.4;
    indices = [1,2];
end
N = length(abc_theta(:,1));
M=11;
abc_theta

%split into bins
edges_vecs = [linspace(prior_params(indices(1))-range/2,prior_params(indices(1))+range/2,M);linspace(prior_params(indices(2))-range/2,prior_params(indices(2))+range/2,M)];
counts = zeros(length(edges_vecs(1,:)));
for k=1:N      
    bin_index = [1,1];  
    for j=1:2
        while edges_vecs(j,bin_index(j))<abc_theta(k,indices(j))
            bin_index(j) = bin_index(j)+1;
        end
    end
    counts(bin_index(1),bin_index(2)) = counts(bin_index(1),bin_index(2))+1;
end
%[Num_in_bins,~] = histc(abc_theta(:,indices),meshgrid(linspace(prior_params(indices(1))-range/2,prior_params(indices(1))+range/2,11),linspace(prior_params(indices(2))-range/2,prior_params(indices(2))+range/2,11)))
posterior = counts/N; 

figure;
contourf(posterior)
figure;
pcolor(posterior)

end