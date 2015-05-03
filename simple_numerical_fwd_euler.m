function simple_numerical_fwd_euler(time_input)
%created JH 30/4/15
%last edit 30/4/15
%solving evolution equation avergaed over theta
%ie just transport equation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L=30; %30 microns
T = time_input; %*60^2; % T in seconds, time_input in hrs
dt = 0.01;
dx = 0.1;
num_t_pts = T/dt+1
num_space_pts = L/dx+2

phi = 58;
speed = 0.4; %microns per second
F1 = (4*(phi/100)-2)/pi; %integral on sin(theta)*T(theta) d theta 
dt/dx*speed*F1

q = zeros(num_space_pts,num_t_pts);

%q(4:10,1) = 1;
q(:,1) = poisspdf(1:num_space_pts,2/dx);

for j=2:num_t_pts

q(2:(num_space_pts-1),j) = 0.5*(q(1:(num_space_pts-2),j-1)+q(1:(num_space_pts-2),j-1)) ...
    + dt/dx*speed*F1/2*(q(1:(num_space_pts-2),j-1) - q(3:(num_space_pts),j-1));
q(num_space_pts,j) = 0;
q(1,j) = q(3,j);

end
close all
figure;
plot(-dx:dx:L,q(:,1),'g-',-dx:dx:L,q(:,num_t_pts))
xlabel('x')
ylabel('probability')
