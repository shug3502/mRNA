function [nu, lambda, phi] = estimators_for_perfect_data(N)

%last edit 15/1/16
%created 15/1/16 
%assume simple model with only one mode of motion, ie AT
%assume we can get perfect knowledge of the data
%That is we observe exactly the positions and times at which the jumps occur

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%set parameters to use for simulating data

    params.input_time = 10;
    params.with_anchoring = 1;
    params.num_modes = 1;
    params.with_plot = 0;
    %Now with data on nurse cells from Alex Davidson
    params.nu1 = 1.16; %speed of RNP complex under active transport [zimyanin e$
    params.nu2 = 0; %ratio between speed for active transport vs diffusion [$
    params.lambda_1=0;   %1/0.13; %transition rate =7.69 [zimyanin et al 2008]
    params.lambda_2 = 0;
    params.omega_1= 0.42;    %1/6*(num_modes>1); %rate of falling off the micro$
    params.omega_2 = 0.84;
    params.phi = 0.58; %percentage of microtubules in posterior direction for b$
    params.x_0=0.5;  %Initially in first compartment, ie. at NPC
    params.Lx = 52; %length of cell in x direction
    params.Ly = 37; %in y direction
    params.Lz = 37; %in z direction
    params.nuc_radius = 10; %radius of nucleus
    params.theta_0 = 0; %initial angle is 0
    params.rc_width = 1;
    params.ellipsoid_boundary = 1;

%Note the rate or switching parameter for one mode is params.omega_1 + params.omega_2

%N=1; %now in imput
paths = zeros(N,10^4,3);
times = zeros(N,10^4);
num_jumps = zeros(N,1);

%generate and store data
for j=1:N
	[~,~, ~, temp_path, temp_times] = velocityjump3D_with_nucleus(params);
	l = numel(temp_times);
	paths(j,1:l,:) = temp_path;
	times(j,1:l) = temp_times;
	num_jumps(j) = l-1; %subtract one for initial point
end

%now create the estimates
%for speed
aux = sqrt(sum(diff(paths,1,2).^2,3))./diff(times,1,2);
aux_nu = zeros(N,1);
for i=1:N	 %have to loop because different numbers of jumps for each sample
	ll = 1;  % num_jumps(i); % can edit so only use 1st jump to get accurate prediction
	aux_nu(i) = sum(aux(i,1:ll),2)./ll;  %estimate for each sample
end
nu = 1/N*sum(aux_nu,1);

%for transition rate
aux = diff(times,1,2);
aux_lambda = zeros(N,1);
for i=1:N
aux_lambda(i) = 1/(sum(aux(i,1:num_jumps(i)),2)./num_jumps(i)); %estimate for each sample
end
lambda = 1/N*sum(aux_lambda,1);
lambda_bounds = lambda*(ones(2,1)+1.96/sqrt(sum(num_jumps))*[-1;1])

%for MT bias
theta = zeros(N,10^4);
aux_phi = zeros(N,1);
d = diff(paths,1,2);
for i=1:N
theta(i,1:num_jumps(i)) = acos(d(i,1:num_jumps(i),1)./sqrt(sum(d(i,1:num_jumps(i),:).^2,3)));  
aux_phi(i) = 1 - sum((theta(i,1:num_jumps(i))>pi/2).*(theta(i,1:num_jumps(i))<3*pi/2),2)./num_jumps(i);
end
phi = 1/N*sum(aux_phi);
