function numerical_solution_of_evolution_eqn_v2(input_time,plot_option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Created 23 April 2015
%Last edit 23 April 2015
%solves evolution equation for 1D velocity jump process model
%for transport of an RNP along a MT by molecular motors
%Assumes: 1 mode of motion, motor can reverse

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
speed = 0.5;
dx = 1;
dv=0.01;
dt = 0.1;
x_init = dx;
v_init = 0.2;
time = 2; %60^2*input_time; %input in hours, time in seconds

phi = 58/100;
lambda = 0;%1/0.13;

x_posterior = 0;
x_anterior = 20; %microns


num_spatial_pts = 1+(x_anterior-x_posterior)/dx %may need to add extra points for the boundaries etc
num_vel_pts = 1+2*speed/dv
num_time_pts = time/dt
p = zeros(num_spatial_pts,num_vel_pts, num_time_pts);
%p_hat = zeros(num_spatial_pts,num_time_pts);
transition_kernel = 2*phi/pi*(asin((-speed:dv:speed)/speed)<pi/2).*(asin((-speed:dv:speed)/speed)>=0) ... 
    + 2*(1-phi)/pi*(asin((-speed:dv:speed)/speed)>=-pi/2).*(asin((-speed:dv:speed)/speed)<0);
transition_kernel = transition_kernel/sum(transition_kernel);
%error('stop NOW!')

%INITIALIZE
p(1+(x_init-x_posterior)/dx,1+round((v_init+speed)/dv),1) = 1; 
t = 1;
while t<(num_time_pts+1)    
    diff_p = repmat(-(-speed:dv:speed),num_spatial_pts-2,1).*(p(3:(num_spatial_pts),1:num_vel_pts,t) ... 
        -p(2:(num_spatial_pts-1),1:num_vel_pts,t))/dx
    p(2:(num_spatial_pts-1),1:num_vel_pts,t+1) = p(2:(num_spatial_pts-1),1:num_vel_pts,t) ...
        + dt*(repmat(-(-speed:dv:speed),num_spatial_pts-2,1).*(p(3:(num_spatial_pts),1:num_vel_pts,t) ... 
        -p(2:(num_spatial_pts-1),1:num_vel_pts,t))/dx - lambda*p(2:(num_spatial_pts-1),1:num_vel_pts,t) ... 
        + lambda*repmat(transition_kernel,num_spatial_pts-2,1).*repmat(sum(p(2:(num_spatial_pts-1),:,t)*dv,2),1,num_vel_pts));
    %BOUNDARY CONDITIONS
    %p(1,:,t+1) = 1;
    p(1,:,t+1) = p(2,:,t+1); %or some other more correct BC
    p(num_spatial_pts,:,t+1) = 0;
    p_hat = sum(p(:,:,t)*dv,2);
    %p(1:5,:,t)
    t=t+1
end

if plot_option
close all
figure;
%p_hat = sum(p(:,:,num_time_pts)*dv,2)
plot(x_posterior:dx:x_anterior,p_hat)
end
toc
end

%   0.369239467973197 for v>=0 
%   0.267380304394384 for v<0
