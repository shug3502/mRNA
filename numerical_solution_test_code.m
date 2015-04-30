function numerical_solution_test_code(input_time,plot_option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Created 24 April 2015
%Last edit 24 April 2015
%solves evolution equation for 1D velocity jump process model
%for transport of an RNP along a MT by molecular motors
%Assumes: 1 mode of motion, motor can reverse
%For debugging

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
speed = 0.5;
dv=0.01;
dt = 0.1;
dx = speed*dt;
x_init = 0.5;
v_init = 0.2;
time = 1; %60^2*input_time; %input in hours, time in seconds

phi = 58/100;
lambda = 0;

x_posterior = 0;
x_anterior = 1; %microns


num_spatial_pts = 1+(x_anterior-x_posterior)/dx %may need to add extra points for the boundaries etc
num_vel_pts = 2
num_time_pts = time/dt
p = zeros(num_spatial_pts,num_vel_pts, num_time_pts);


transition_kernel = 1;

%INITIALIZE
p(1+(x_init-x_posterior)/dx,1,1) = 1; 
t = 1;
while t<(num_time_pts+1)    
    
    p(2:(num_spatial_pts-1),1,t+1) = p(1:(num_spatial_pts-2),1,t) ...
        + dt*( - lambda*p(2:(num_spatial_pts-1),1,t) ... 
        + lambda*transition_kernel.*p(2:(num_spatial_pts-1),2,t));
    p(2:(num_spatial_pts-1),2,t+1) = p(3:(num_spatial_pts),2,t) ...
        + dt*(- lambda*p(2:(num_spatial_pts-1),2,t) ... 
        + lambda*transition_kernel.*p(2:(num_spatial_pts-1),1,t));
    
%     p(2:(num_spatial_pts-1),1,t+1) = p(2:(num_spatial_pts-1),1,t) ...
%         + dt*(-speed.*(p(3:(num_spatial_pts),1,t) ... 
%         -p(1:(num_spatial_pts-2),1,t))/(2*dx) - lambda*p(2:(num_spatial_pts-1),1,t) ... 
%         + lambda*transition_kernel.*p(2:(num_spatial_pts-1),2,t));
%     p(2:(num_spatial_pts-1),2,t+1) = p(2:(num_spatial_pts-1),2,t) ...
%         + dt*(speed.*(p(3:(num_spatial_pts),2,t) ... 
%         -p(1:(num_spatial_pts-2),2,t))/(2*dx) - lambda*p(2:(num_spatial_pts-1),2,t) ... 
%         + lambda*transition_kernel.*p(2:(num_spatial_pts-1),1,t));
%     

%  p(2:(num_spatial_pts-1),1,t+1) = p(2:(num_spatial_pts-1),1,t) ...
%         + dt*(speed.*(p(3:(num_spatial_pts),1,t) - 2*p(2:(num_spatial_pts-1),1,t) ... 
%         +p(1:(num_spatial_pts-2),1,t))/(dx^2) - lambda*p(2:(num_spatial_pts-1),1,t) ... 
%         + lambda*transition_kernel.*p(2:(num_spatial_pts-1),2,t));
%     p(2:(num_spatial_pts-1),2,t+1) = p(2:(num_spatial_pts-1),2,t) ...
%         + dt*(speed.*(p(3:(num_spatial_pts),2,t) - 2*p(2:(num_spatial_pts-1),2,t) ... 
%         +p(1:(num_spatial_pts-2),2,t))/(dx^2) - lambda*p(2:(num_spatial_pts-1),2,t) ... 
%         + lambda*transition_kernel.*p(2:(num_spatial_pts-1),1,t));

    %BOUNDARY CONDITIONS
    p(1,:,t+1) = p(2,:,t+1); %or some other more correct BC
    p(num_spatial_pts,:,t+1) = 0;
    size(p)
    p_hat = sum(p(:,:,t),2)
    t=t+1
end

if plot_option
close all
figure;
p_hat = sum(p(:,:,num_time_pts),2)
length(x_posterior:dx:x_anterior)
length(p_hat)
plot(x_posterior:dx:x_anterior,p_hat)
end
toc
end
