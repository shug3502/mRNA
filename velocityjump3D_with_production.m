function [q_raw,q_distn] = velocityjump3D_with_production(params)
%Created 18/3/16
%Last edit 4/5/16
%Based on 2D version
%Uses the new boundary conditions
%3D velocity jump model for transport of an RNP along a MT by molecular motors
%Includes production in nucleus
%Assumes: 1 mode of motion, diffusion is neglected
%Uses biologically motivated parameters
%Edited to transitions between states are clearer and nucleus is included
%Edited so that parameter input is much nicer

%with_anchoring is logical to include anhoring at site of localization or reflection there
%with_plot is logical to include plotting or not
%my_seed is the random seed to be used
%num_modes is for the number of modes of motion ie only AT or with
%diffusion also
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin ~= 1
    %fprintf('default parameters used\n')
    params.input_time = 1;
    params.with_anchoring = 1;
    params.num_modes = 1;
    params.with_plot = 0;
    %Now with data on nurse cells from Alex Davidson
    params.nu1 = 1.16; %speed of RNP complex under active transport [zimyanin et al 2008]
    params.nu2 = 0.80; %ratio between speed for active transport vs diffusion [zimyanin et al 2008]
    params.lambda_1=0;   %1/0.13; %transition rate =7.69 [zimyanin et al 2008]
    params.lambda_2 = 0.11;
    params.omega_1= 0.42;    %1/6*(num_modes>1); %rate of falling off the microtubule [zimyanin et al 2008] since average track length 2.4 - 2.8 microns -> average jump for 6s -> rate 1/6
    params.omega_2 = 0.84;
    params.phi = 0.58; %percentage of microtubules in posterior direction for biased angle distn [parton et al 2011]
    params.x_0=0.5;  %Initially in first compartment, ie. at NPC
    params.Lx = 52; %length of cell in x direction
    params.Ly = 37; %in y direction
    params.Lz = 37; %in z direction
    params.nuc_radius = 10; %radius of nucleus
    params.theta_0 = 0; %initial angle is 0
    params.rc_width = 1;
    params.ellipsoid_boundary = 1;
    params.gamma=0.1;
end

tic
t_end = params.input_time;
N=ceil(4*t_end*60^2*params.gamma); %guess for number of particles (guess larger than likely value)

l_n = 5000; %storage for the jumps and paths
delx=1;
L=params.Lx;
time_vec = (0:0.5:t_end);
t=time_vec*60^2;
l_t = length(time_vec);
jumps = zeros(N,l_n);
q_distn = zeros(L/delx+1,l_t);
q_raw = q_distn;
xpos = zeros(N,l_n);
xpos_discrete_time = -1*ones(N,l_n);

time = 0;
num_particles=0;
while time<t_end*60^2    
    num_particles=num_particles+1; %count num particles produced   
    params.input_time = t_end-time/60^2;
    [is_anchored, ~,~, path, jump_times] = velocityjump3D_with_nucleus(params);
    xpos_temp = path(:,1);
    jump_temp = jump_times + time;
    jumps(num_particles,1:length(jump_temp)) = jump_temp;
    xpos(num_particles,1:length(xpos_temp)) = xpos_temp;
    for w=1:l_t
        find_index = find(t(w)<jumps(num_particles,:),1,'first');
        if isempty(xpos(num_particles,find_index)) && is_anchored
            xpos_discrete_time(num_particles,w) = L;
        else
            if find_index>1
            xpos_discrete_time(num_particles,w) =  min(xpos(num_particles,find_index-1),L);
            end
        end
    end
    
    %update time after one reaction (production event)
    tor = -1/params.gamma*log(rand(1));
    time = time + tor;
end

xpos_discrete_time = xpos_discrete_time(1:num_particles,:);
for w=1:l_t
    %split into bins
    [Num_in_bins,~] = histc(xpos_discrete_time(:,w),0:delx:L);
%     if sum(Num_in_bins) ~= num_particles
%         xpos_discrete_time(:,w)
%         ename = sprintf('ABC_errors.txt');
%         fileIDerror = fopen(ename,'w');
%         fprintf(fileIDerror,'%f \n',xpos_discrete_time);
%         fclose('all');
%         error('outside of 0:L');
%     end
    q_distn(:,w) = Num_in_bins/sum(Num_in_bins); %/num_particles/delx; %estimate of q at time T
    q_raw(:,w) = Num_in_bins;
end
if params.with_plot
%plot with a log scale
figure;
subplot(2,1,1)
imagesc(time_vec(1:9), 0:L, flipud(log(q_raw(:,1:9)))) %plot initial period
colormap(gray)
set(gca, 'fontsize', 20);
yticklabels = 0:10:52;
yticks = linspace(1, size(q_raw, 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)))
xlabel('Time')
ylabel('Position')
subplot(2,1,2)
imagesc(time_vec, 0:L, flipud(log(q_raw))) % plot whole time course
colormap(gray)
set(gca, 'fontsize', 20);
yticklabels = 0:10:52;
yticks = linspace(1, size(q_raw, 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)))
xlabel('Time')
ylabel('Position')

%plot on a normal scale 
figure;
subplot(2,1,1)
imagesc(time_vec(1:9), 0:L, flipud(q_raw(:,1:9))) % plot initial period
colormap(gray)
set(gca, 'fontsize', 20);
yticklabels = 0:10:52;
yticks = linspace(1, size(q_raw, 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)))
xlabel('Time')
ylabel('Position')
subplot(2,1,2)
imagesc(time_vec, 0:L, flipud(q_raw))
colormap(gray)
set(gca, 'fontsize', 20);
yticklabels = 0:10:52;
yticks = linspace(1, size(q_raw, 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)))
xlabel('Time')
ylabel('Position')
end
toc
end

