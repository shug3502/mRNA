function posterior = bin_posterior(abc_theta,prior_params,range)
%prior params should be a vector with the centres of the uniform
%distribution for the prior
%range is the spread of the normal prior. Assumed this is the same for all
%the parameters

if nargin~=3
fname = sprintf('Simplest_ABC_with_moments_output2.txt');
fileID = fopen(fname,'r');
abc_theta = fscanf(fileID, '%f');
N = length(abc_theta)/3;
abc_theta = reshape(abc_theta,N,3);
fclose('all');
prior_params = [1.16,0.42,0.58];
range = 0.4;
end
N = length(abc_theta(:,1));
for j=1:3
    %split into bins
    [Num_in_bins,~] = histc(abc_theta(:,j),linspace(prior_params(j)-range/2,prior_params(j)+range/2,11));
    density_estimate(:,j) = Num_in_bins/N; %estimate of q at time T
end
posterior = density_estimate;

figure;
contourf(posterior)

end