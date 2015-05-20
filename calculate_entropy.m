function entropy = calculate_entropy(abc_theta)
%calculate entropy via k-l divergence
%High entropy means low information

if nargin ~= 1
    fprintf('default parameters taken from file\n')
fname = sprintf('Simplest_ABC_with_moments_output.txt');
fileID = fopen(fname,'r');
abc_theta = fscanf(fileID, '%f');
fclose('all');
abc_theta = reshape(abc_theta,length(abc_theta)/3,3);
end

bin_width = 0.01;
param_bins = [0.7:bin_width:1.1; %nu_1
  7.5:bin_width:7.9; %lambda
  0.38:bin_width:0.78]; %phi
distn_theta = zeros(size(param_bins));
for j=1:3
 [distn_theta(j,:),~] = histc(abc_theta(:,j),param_bins(j,:));
end

entropy = 0;
for i=1:3
entropy = entropy + kldiv(param_bins(1,:)+bin_width/2, distn_theta(i,:)/sum(distn_theta(i,:))+eps,ones(size(distn_theta(i,:)))/length(distn_theta));
end

end
