function numerical_solution_of_evolution_eqn(input_time,plot_option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Created 23 April 2015
%Last edit 23 April 2015
%solves evolution equation for 1D velocity jump process model
%for transport of an RNP along a MT by molecular motors
%Assumes: 1 mode of motion, motor can reverse

%Too slow!!!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
speed = 0.4;
dx = 1;
dv=0.01;
dt = 1;
x_init = dx;
v_init = 0;
time = 60^2*input_time; %input in hours, time in seconds

phi = 58;
lambda = 1/0.13;

x_posterior = 0;
x_anterior = 80; %microns


num_spatial_pts = 1+(x_anterior-x_posterior)/dx %may need to add extra points for the boundaries etc
num_vel_pts = 1+2*speed/dv
num_time_pts = time/dt
p = zeros(num_spatial_pts,num_vel_pts, num_time_pts);
%p_hat = zeros(num_spatial_pts,num_time_pts);

%INITIALIZE
p(1+(x_init-x_posterior)/dx,1+(v_init+speed)/dv,1) = 1; 
t = 1;
while t<num_spatial_pts
    
    for i=2:num_spatial_pts-1
        p_hat = sum(p*dv,2);
        for j=1:num_vel_pts
            p(i,j,t+1) = p(i,j,t) + dt*(-((j-1)*dv-speed)*(p(i+1,j,t)-p(i,j,t))/dx - lambda*p(i,j,t) + lambda*transition_kernel((j-1)*dv-speed,speed,phi)*p_hat(i,t));
        end
    end
    p(1,:,t+1) = p(2,:,t); %or some other more correct BC
    p(num_spatial_pts,:,t+1) = 0;
    t=t+1
end

if plot_option
close all
figure;
p_hat = sum(p(:,:,num_time_pts)*dv,2);
plot(x_anterior:dx:x_posterior,p_hat)
end
toc
end

function T = transition_kernel(v_new,speed,phi)
T = 2*phi/pi*(asin(v_new/speed)<pi/2)*(asin(v_new/speed)>0) + 2*(1-phi)/pi*(asin(v_new/speed)>-pi/2)*(asin(v_new/speed)<0); 

end
