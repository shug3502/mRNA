function mfpt = VJP_repeated_fwd_simulate(theta,hyperparams,repeats)
%JH 08/10/18
%wrap velocity jump process 3d for ABC with mfpt as summary statistic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hyperparams.nu1 = theta(1); %speed of RNP complex under active transport [zimyanin et al 2008]
hyperparams.nu2 = theta(2); %ratio between speed for active transport vs diffusion [zimyanin et al 2008]
hyperparams.lambda_1 = 0;   %transition rate [zimyanin et al 2008]
hyperparams.lambda_2 = theta(3);
hyperparams.omega_1= theta(4);    %rate of falling off the microtubule [zimyanin et al 2$
hyperparams.omega_2 = theta(5);
hyperparams.phi = theta(6); %percentage of microtubules in posterior direction for biased angle distn [parto$

anchor_times = zeros(repeats,1);
for j=1:repeats
    [~, anchor_times(j), ~, ~, ~] = velocityjump3D_with_nucleus_and_new_BCs(hyperparams);
end
mfpt = mean(anchor_times);

end
