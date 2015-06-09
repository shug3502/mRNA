function [entropy, quality, proportion_in_interval] = compare_posterior_qualities(N,num_runs)
%compare quality of posterior from ABC with ABC_APMC and ABC_weights
%dependent on distance

prior_params = [1.16, 0.8, 0.11, 0.42, 0.84, 0.58, 0.5, 0];
real_params = prior_params; %if real params used are different then MUST change this!
real_params(1) = 1;
range = 0.8;
indices = [1,4,6];
alpha = 0.5;
% abc_theta = zeros(N*num_runs*alpha,3);
% abc_weights = zeros(N*num_runs*alpha,1);
%[posterior,quality_of_posterior] = multi_dim_bin_posterior(abc_theta,abc_weights,prior_params,real_params,range,indices,1);

%num_runs = 20;



for k=[0,1] %0 is with spatial distribution; 1 is mfpt etc.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    proportion_in_interval = 0;
    entropy = 0;
    quality = 0;
    abc_theta = [];
    abc_weights = [];
    
    for b=1:num_runs
        [abc_theta_temp,abc_weights_temp,entropy_temp,quality_temp,contained_in_pred_interval] =ABC_RS_knn(N,k,b);
        abc_theta = [abc_theta; abc_theta_temp];
        abc_weights = [abc_weights; abc_weights_temp];       
%         proportion_in_interval = proportion_in_interval+contained_in_pred_interval;
%         entropy = entropy+entropy_temp;
%         quality = quality+quality_temp;
    end
    
%     proportion_in_interval = proportion_in_interval/num_runs;
%     entropy = entropy/num_runs;
%     quality = quality/num_runs;
    [entropy, quality, proportion_in_interval]= post_processing(abc_theta,abc_weights,prior_params,real_params,repmat(range,1,3),indices);
    
    fname = sprintf('Comparison_of_ABC_methods_N100_summarys_edited.txt');
    fileID = fopen(fname,'a');
    fprintf(fileID,'%d %f %f %f %f \n',k,entropy,quality(1),quality(2),proportion_in_interval);
    fclose('all');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    proportion_in_interval = 0;
    entropy = 0;
    quality = 0;
    abc_theta = [];
    abc_weights = [];
    
    for b=1:num_runs
        [abc_theta,abc_weights,entropy_temp,quality_temp,contained_in_pred_interval] =ABC_Weights_depend_on_dist(N,k,b);
        abc_theta = [abc_theta; abc_theta_temp];
        abc_weights = [abc_weights; abc_weights_temp];       
%         proportion_in_interval = proportion_in_interval+contained_in_pred_interval;
%         entropy = entropy+entropy_temp;
%         quality = quality+quality_temp;
    end
    
%     proportion_in_interval = proportion_in_interval/num_runs;
%     entropy = entropy/num_runs;
%     quality = quality/num_runs;
    [entropy, quality, proportion_in_interval]= post_processing(abc_theta,abc_weights,prior_params,real_params,repmat(range,1,3),indices);
    
    fileID = fopen(fname,'a');
    fprintf(fileID,'%d %f %f %f %f \n',k,entropy,quality(1),quality(2),proportion_in_interval);
    fclose('all');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    %Now for APMC
    proportion_in_interval = 0;
    entropy = 0;
    quality = 0;
    abc_theta = [];
    abc_weights = [];
    
    for b=1:num_runs
        [abc_theta,abc_weights,entropy_temp,quality_temp,contained_in_pred_interval] =ABC_APMC(N,k,b);
        abc_theta = [abc_theta; abc_theta_temp];
        abc_weights = [abc_weights; abc_weights_temp];       
%         proportion_in_interval = proportion_in_interval+contained_in_pred_interval;
%         entropy = entropy+entropy_temp;
%         quality = quality+quality_temp;
    end
    
%     proportion_in_interval = proportion_in_interval/num_runs;
%     entropy = entropy/num_runs;
%     quality = quality/num_runs;
    [entropy, quality, proportion_in_interval]= post_processing(abc_theta,abc_weights,prior_params,real_params,repmat(range,1,3),indices);
    
    fileID = fopen(fname,'a');
    fprintf(fileID,'%d %f %f %f %f \n',k,entropy,quality(1),quality(2),proportion_in_interval);
    fclose('all');
    
end

end
