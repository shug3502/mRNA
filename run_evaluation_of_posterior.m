function [entropy, quality, proportion_in_interval] = run_evaluation_of_posterior(N,num_runs)
%evaluate posterior for single method

for k=[0,1] %0 is with spatial distribution; 1 is mfpt etc.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
proportion_in_interval = 0;
entropy = 0;
quality = 0;

for b=1:num_runs
    [~,~,entropy_temp,quality_temp,contained_in_pred_interval] =ABC_RS_knn(N,k,b);
%     abc_theta((1+(b-1)*N*alpha):(N*b*alpha),:) = abc_theta_temp;
%     abc_weights((1+(b-1)*N*alpha):(N*b*alpha),:) = abc_weights_temp;    
    proportion_in_interval = proportion_in_interval+contained_in_pred_interval;
    entropy = entropy+entropy_temp;
    quality = quality+quality_temp;
end

proportion_in_interval = proportion_in_interval/num_runs;
entropy = entropy/num_runs;
quality = quality/num_runs;


fname = sprintf('Comparison_of_ABC_methods_prior.txt');
fileID = fopen(fname,'a');
fprintf(fileID,'%d %f %f %f %f \n',k,entropy,quality(1),quality(2),proportion_in_interval);
fclose('all');
end

end