function backward_euler_advection

clear;
%problem parameters
L = 1.0;
Tmax = 1.0;
c = 1.0;

%parameters for the solution
maxt = 350; %number of time steps
dt = Tmax/maxt;
n = 300; %number of space steps
nint = 50; %the wave front
dx = L/n;
b = c*dt/(2.*dx);

%initial value of the function u (amplitude of the wave)
for i = 1:(n+1)
    if i < nint
        u(i,1) = 1.;
    else
        u(i,1)=0.;
    end
    x(i) = (i-1)*dx;
end

%Value of amplitude at boundary
for k=1:maxt+1
    u(1,k) = 1.;
    u(n+1,k) = 0.;
    time(k) = (k-1)*dt;
end

%implement the explicit method
for k=1:maxt
    for i=2:n
        u(i,k+1) = 0.5*(u(i+1,k)+u(i-1,k)) - b*(u(i+1,k)-u(i-1,k));
    end
end

plot(x,u(:,1),'-',x,u(:,10),'-',x,u(:,50),'-',x,u(:,100),'-');

