function [entropy, quality, proportion_in_interval] = compare_posterior_qualities_single_version(N,num_runs,which_method)
%compare quality of posterior from ABC with ABC_APMC and ABC_weights
%dependent on distance

[~, real_params, ~, ~, ...
    ~, prior_params, indices, sigma, ~] = initialise_parameters(0,0);


% prior_params = [5, 0.8, 0.11, 0.42, 0.84, 0.58, 0.5, 0];
% real_params = prior_params; %if real params used are different then MUST change this!
% real_params(1) = 1.16;
% sigma = 10;
% indices = 1; %[1,4,6];
% alpha = 0.5;
% % abc_theta = zeros(N*num_runs*alpha,3);
% % abc_weights = zeros(N*num_runs*alpha,1);
% %[posterior,quality_of_posterior] = multi_dim_bin_posterior(abc_theta,abc_weights,prior_params,real_params,range,indices,1);

%num_runs = 20;

fname = sprintf('Comparison_of_ABC_methods_N1000_RS.txt');
gname = sprintf('ABC_theta_store_N1000_RS.txt');

for k=[0, 1] %0 is with spatial distribution; 1 is mfpt etc.
   
    if which_method==1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    proportion_in_interval = 0;
    entropy = 0;
    quality = 0;
    abc_theta = [];
    abc_weights = [];
    
    for b=1:num_runs
        [abc_theta_temp,abc_weights_temp,entropy_temp,quality_temp,contained_in_pred_interval] =ABC_RS_knn(N,k,b,0);
        abc_theta = [abc_theta; abc_theta_temp];
        abc_weights = [abc_weights; abc_weights_temp];       
%         proportion_in_interval = proportion_in_interval+contained_in_pred_interval;
%         entropy = entropy+entropy_temp;
%         quality = quality+quality_temp;
    end
    
%     proportion_in_interval = proportion_in_interval/num_runs;
%     entropy = entropy/num_runs;
%     quality = quality/num_runs;
    [entropy, quality, proportion_in_interval]= post_processing(abc_theta,abc_weights,prior_params,real_params,sigma,indices);
    
     
    fileID = fopen(fname,'a');
    fprintf(fileID,'%d %f %f %f %f \n',k,entropy,quality(1),quality(2),proportion_in_interval);
    
    fileID2 = fopen(gname,'a');
    fprintf(fileID2,'%d \n \n %f \n',k,abc_theta);
    fclose('all');
    
    elseif which_method == 2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    proportion_in_interval = 0;
    entropy = 0;
    quality = 0;
    abc_theta = [];
    abc_weights = [];
    
    for b=1:num_runs
        [abc_theta_temp,abc_weights_temp,entropy_temp,quality_temp,contained_in_pred_interval] =ABC_Weights_depend_on_dist(N,k,b,0);
        abc_theta = [abc_theta; abc_theta_temp];
        abc_weights = [abc_weights; abc_weights_temp];       
%         proportion_in_interval = proportion_in_interval+contained_in_pred_interval;
%         entropy = entropy+entropy_temp;
%         quality = quality+quality_temp;
    end
    
%     proportion_in_interval = proportion_in_interval/num_runs;
%     entropy = entropy/num_runs;
%     quality = quality/num_runs;
    [entropy, quality, proportion_in_interval]= post_processing(abc_theta,abc_weights,prior_params,real_params,sigma,indices);
    
    fileID = fopen(fname,'a');
    fprintf(fileID,'%d %f %f %f %f \n',k,entropy,quality(1),quality(2),proportion_in_interval);
    fileID2 = fopen(gname,'a');
    fprintf(fileID2,'%d %f \n',k,abc_theta);
    fclose('all');
    
    elseif which_method == 3
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %Now for APMC
    proportion_in_interval = 0;
    entropy = 0;
    quality = 0;
    abc_theta = [];
    abc_weights = [];
    
    for b=1:num_runs
        [abc_theta_temp,abc_weights_temp,entropy_temp,quality_temp,contained_in_pred_interval] =ABC_APMC(N,k,b,0);
        abc_theta = [abc_theta; abc_theta_temp];
        abc_weights = [abc_weights; abc_weights_temp];       
%         proportion_in_interval = proportion_in_interval+contained_in_pred_interval;
%         entropy = entropy+entropy_temp;
%         quality = quality+quality_temp;
    end
    
%     proportion_in_interval = proportion_in_interval/num_runs;
%     entropy = entropy/num_runs;
%     quality = quality/num_runs;
    [entropy, quality, proportion_in_interval]= post_processing(abc_theta,abc_weights,prior_params,real_params,sigma,indices);
    
    fileID = fopen(fname,'a');
    fprintf(fileID,'%d %f %f %f %f \n',k,entropy,quality(1),quality(2),proportion_in_interval);
    fileID2 = fopen(gname,'a');
    fprintf(fileID2,'%d %f \n',k,abc_theta);
    fclose('all');
    end
end

end
