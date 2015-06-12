function entropy = calculate_entropy(abc_theta,prior_params,prior_sigma,p_indices)
%calculate entropy via k-l divergence
%High entropy means low information

if nargin ~= 4
    fprintf('default parameters taken from file\n')
fname = sprintf('Simplest_ABC_with_moments_output.txt');
fileID = fopen(fname,'r');
abc_theta = fscanf(fileID, '%f');
fclose('all');
abc_theta = reshape(abc_theta,length(abc_theta)/3,3);
end

%bin_width = 0.01;
M=21;
%param_bins = zeros(3,max(prior_sigma)/bin_width+1);
%(prior_params(p_indices(i))-prior_sigma(i)/2):bin_width:(prior_params(p_indices(i))+prior_sigma(i)/2); 
for i=1:length(p_indices)
param_bins(i,:) = linspace((prior_params(p_indices(i))-prior_sigma(i)/2),(prior_params(p_indices(i))+prior_sigma(i)/2),M); 
end

distn_theta = zeros(size(param_bins));
for j=1:length(p_indices)
 [distn_theta(j,:),~] = histc(abc_theta(:,j),param_bins(j,:));
end
entropy = 0;
for i=1:length(p_indices)
entropy = entropy + kldiv(param_bins(i,:)+prior_sigma(i)/(2*M), distn_theta(i,:)/sum(distn_theta(i,:))+eps,ones(size(distn_theta(i,:)))/length(distn_theta));
end

end
