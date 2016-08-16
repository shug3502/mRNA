function output = mRNA_repeats_wrapper(coords,theta,params)
%wraps mRNA_3D_GFP function, and allows multiple repeats
%intended for mixture model fitting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<1
    coords = get_coords('thresholded');
    theta = [1.16,0.8,0.58,0.42,0.84,0.11];
    t_end = 1; %3*4.5*10^-3;
    dt=0.1; %1.5*10^-4;
    params.with_anchoring = 1;
    params.num_modes = 2;
    params.nu1 = theta(1); %speed of RNP complex under active transport [zimyanin et al 2008]
    params.nu2 = theta(2); %ratio between speed for active transport vs diffusion [zimyanin et al 2008]
    params.lambda_1=0;   %1/0.13; %transition rate =7.69 [zimyanin et al 2008]
    params.lambda_2 = theta(6);
    params.omega_1= theta(4);    %1/6*(num_modes>1); %rate of falling off the microtubule [zimyanin et al 2008] since average track length 2.4 - 2.8 microns -> average jump for 6s -> rate 1/6
    params.omega_2 = theta(5);
    params.phi = theta(3); %percentage of microtubules in posterior direction for biased angle distn [parton et al 2011]
    params.Lx = 52; %length of cell in x direction
    params.Ly = 37; %in y direction
    params.Lz = 37; %in z direction
    params.nuc_radius = 10; %radius of nucleus
    params.theta_0 = 0; %initial angle is 0
    params.rc_width = 1;
    params.ellipsoid_boundary = 1;
    params.sample_times = 0:dt:t_end; %roughly the same as time series for thresholded image
    params.num_hist_bins=21;
    params.zslice = [-1,1];
    params.path_summary_stat=1;
    params.repeats = 10;
end
sz = params.num_hist_bins*numel(params.sample_times);
b = zeros(1,sz,params.repeats);
for k=1:params.repeats
    b(1,:,k) = reshape(mRNA_3D_GFP(coords,params),[],1);
end
output = b;
