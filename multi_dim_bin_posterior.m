function [posterior,quality_of_posterior] = multi_dim_bin_posterior(abc_theta,abc_weights,prior_params,real_params,range,indices,plot_option)
%prior params should be a vector with the centres of the uniform
%distribution for the prior
%range is the spread of the normal prior. Assumed this is the same for all
%the parameters
%computes multidimensional histogram

if nargin~=7
    fprintf('Using defaults \n');
    fname = sprintf('ABC_theta_read_from_R.dat');   %('Simplest_ABC_with_moments_output2.txt');
%     fileID = fopen(fname,'r');
%     abc_theta = fscanf(fileID, '%f')
%     N = length(abc_theta)/3
    %abc_theta = reshape(abc_theta,N,3);
%     fclose('all');
    abc_theta = load(fname);
    abc_weights = ones(length(abc_theta),1);
    
    [~, real_params,~, ~, ~, prior_params, indices, range, ~] = initialise_parameters(0,0);
    
%     prior_params = [5, 5, 5, 1.5, 1.5, 0.5, 0.5, 0];
%     real_params = prior_params;
%     range =  [10, 10, 10, 3, 3, 1]; 
%     indices = 1:6; %[1,4,6];
    plot_option = 1;
end
N = length(abc_theta(:,1));
M=21;
posterior = 1;

%split into bins
edges_vecs = zeros(length(indices),M);
for i=1:length(indices)
    edges_vecs(i,1:M) = linspace(prior_params(indices(i))-range(indices(i))/2,prior_params(indices(i))+range(indices(i))/2,M);
end
% [count, ~, mid, ~] = histcn(abc_theta, edges_vecs(1,:),edges_vecs(2,:),edges_vecs(3,:),'AccumData',abc_weights);
% posterior = count/N;


if plot_option
    var_strings = {'\nu_1','\nu_2','\lambda','\omega_1','\omega_2','\phi'};
    posterior_mode = zeros(1,length(indices));
    figure;
    for j=1:length(indices)
        p_string = sprintf('param %d',j);
        subplot(ceil(length(indices)/2),2,j)
        bincounts = histc(abc_theta(:,j),edges_vecs(j,:));
        [~,mode_index] = max(bincounts);
        posterior_mode(j) = (edges_vecs(j,mode_index)+edges_vecs(j,mode_index+1))/2;
        bar(edges_vecs(j,:),bincounts); %,'histc');
        hold on
        line([real_params(indices(j)),real_params(indices(j))],[0,max(bincounts)],'LineWidth',3,'Color','g');
        xlabel(var_strings(j));
        
    end
    print('Auto_hist_plot','-dpng');
    posterior_mode
    posterior_mean = mean(abc_theta,1);
    quality_of_posterior = [sqrt(sum((posterior_mode - real_params(indices)).^2)),sqrt(sum((posterior_mean - real_params(indices)).^2))]; %note this currently does not depend equally on all params due to scaling of each param
    
    figure;
    fig_index=0;
    for i=1:length(indices)
        for j=1:length(indices)
            if i<j
                fig_index = fig_index+1;
                subplot(ceil(nchoosek(length(indices),2)/3),3,fig_index);
                p_string = sprintf('param %d',i);
                q_string = sprintf('param %d',j);
                [count, ~, mid, ~] = histcn(abc_theta(:,[i,j]), edges_vecs(i,:),edges_vecs(j,:),'AccumData',abc_weights);
                imagesc(mid{1:2},count(:,:,ceil(end/4))')
                xlabel(var_strings(i))
                ylabel(var_strings(j))
            end
        end
    end
    print('Auto_heatmap_plot','-dpng');
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