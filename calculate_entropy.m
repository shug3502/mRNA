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

bin_width = 0.01;
%param_bins = zeros(3,max(prior_sigma)/bin_width+1);
for i=1:3
param_bins(i,:) = (prior_params(p_indices(i))-prior_sigma(i)/2):bin_width:(prior_params(p_indices(i))+prior_sigma(i)/2); 
end

distn_theta = zeros(size(param_bins));
for j=1:3
 [distn_theta(j,:),~] = histc(abc_theta(:,j),param_bins(j,:));
end
entropy = 0;
for i=1:3
entropy = entropy + kldiv(param_bins(i,:)+bin_width/2, distn_theta(i,:)/sum(distn_theta(i,:))+eps,ones(size(distn_theta(i,:)))/length(distn_theta));
end

end
