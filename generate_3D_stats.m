%generate stats for 3D VJump model
N=100;

%default parameters
params.input_time = 40;
params.with_anchoring = 1;
params.num_modes = 2;
params.with_plot = 0;
params.nu1 = 1.16; %speed of RNP complex under active transport [zimyanin et al 2008]
params.nu2 = 0.80; %ratio between speed for active transport vs diffusion [zimyanin et al 2008]
params.lambda_1=0;   %1/0.13; %transition rate =7.69 [zimyanin et al 2008]
params.lambda_2 = 0.11;
params.omega_1= 0.42;    %1/6*(num_modes>1); %rate of falling off the microtubule [zimyanin et al 2008] since average track length 2.4 - 2.8 microns -> average jump for 6s -> rate 1/6
params.omega_2 = 0.84;
params.phi = 0.58; %percentage of microtubules in posterior direction for biased angle distn [parton et al 2011]
params.Lx = 52; %length of cell in x direction
params.Ly = 37; %in y direction
params.Lz = 37; %in z direction
params.nuc_radius = 10; %radius of nucleus
params.theta_0 = 0; %initial angle is 0
params.rc_width = 1;

%with a cuboid boundary
params.ellipsoid_boundary = 0;

tic; cube_anch_times = zeros(N,1); 
for j=1:N 
    [is_anchored, anchoring_time, final_position, path, jump_times] = velocityjump3D_with_nucleus(params); 
    cube_anch_times(j) = anchoring_time/60^2; 
	anchoring_time - jump_times(end)
end 
cube_anch_times
toc
fprintf('For cuboid boundary we have: mean %f, median %f ,sd %f \n',mean(cube_anch_times),median(cube_anch_times),std(cube_anch_times));
figure; 
hist(cube_anch_times,20);
title('MFPT distribution with cuboid boundary');
xlabel('MFPT');
ylabel('Frequency');

% now with an ellipsoid boundary
params.ellipsoid_boundary=1;

s = (3/(pi*4))^(1/3); %scaling factor so that volumes of domains are equal
params.Lx = 54*s;
params.Ly = 37*s;
params.Lz = 37*s;

tic; ellipse_anch_times = zeros(N,1); 
for j=1:N
    [is_anchored, anchoring_time, final_position, path, jump_times] = velocityjump3D_with_nucleus(params); 
    ellipse_anch_times(j) = anchoring_time/60^2; 
end
toc
fprintf('For ellipsoid boundary we have: mean %f, median %f, sd %f \n',mean(ellipse_anch_times),median(ellipse_anch_times), std(ellipse_anch_times));
figure; 
hist(ellipse_anch_times,20);
title('MFPT distribution with ellipsoid boundary');
xlabel('MFPT');
ylabel('Frequency');
