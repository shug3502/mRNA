function posterior = multi_dim_bin_posterior(abc_theta,abc_weights,prior_params,real_params,range,indices,plot_option)
%prior params should be a vector with the centres of the uniform
%distribution for the prior
%range is the spread of the normal prior. Assumed this is the same for all
%the parameters
%computes multidimensional histogram

if nargin~=7
    fprintf('Using defaults \n');
    fname = sprintf('ABC_APMC_output2.txt');   %('Simplest_ABC_with_moments_output2.txt');
    fileID = fopen(fname,'r');
    abc_theta = fscanf(fileID, '%f');
    N = length(abc_theta)/3;
    abc_theta = reshape(abc_theta,N,3);
    fclose('all');
    abc_weights = ones(length(abc_theta),1);
    abc_weights(1)=1;
    prior_params = [1.16, 0.8, 0.11, 0.42, 0.84, 0.58, 0.5, 0];
    real_params = prior_params;
    range = 0.4;
    indices = [1,4,6];
    plot_option = 1;
end
N = length(abc_theta(:,1));
M=21;

%split into bins
edges_vecs = [linspace(prior_params(indices(1))-range/2,prior_params(indices(1))+range/2,M);
    linspace(prior_params(indices(2))-range/2,prior_params(indices(2))+range/2,M);
    linspace(prior_params(indices(3))-range/2,prior_params(indices(3))+range/2,M)];

[count, ~, mid, ~] = histcn(abc_theta, edges_vecs(1,:),edges_vecs(2,:),edges_vecs(3,:),'AccumData',abc_weights);
posterior = count/N;

if plot_option
    figure;
    for j=1:3
        p_string = sprintf('param %d',j);
        subplot(3,1,j)
        bincounts = histc(abc_theta(:,j),edges_vecs(j,:));
        bar(edges_vecs(j,:),bincounts,'histc');
        hold on
        line([real_params(indices(j)),real_params(indices(j))],[0,max(bincounts)],'LineWidth',3,'Color','g');
        xlabel(p_string);
    end
%     figure;
%     imagesc(mid{1:2},posterior(:,:,ceil(end/2)));
%     
%     figure;
%     pcolor(0.5*(edges_vecs(1,1:M-1)+edges_vecs(1,2:M)),(edges_vecs(2,1:M-1)+0*edges_vecs(2,1:M-1)), mean(posterior,3));
%     
%     figure;
%     contourf(0.5*(edges_vecs(1,1:M-1)+edges_vecs(1,2:M)),(edges_vecs(2,1:M-1)+0*edges_vecs(2,1:M-1)), mean(posterior,3));
end
end