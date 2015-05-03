function pdex1
close all

m = 0;
x = linspace(0,1,50);
t = linspace(0,2,50);

options = odeset('RelTol',1e-5, 'AbsTol',1e-5);
sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t,options);
% Extract the first solution component as u.
u = sol(:,:,1);

% A surface plot is often a good way to study a solution.
surf(x,t,u) 
title('Numerical solution computed with 20 mesh points.')
xlabel('Distance x')
ylabel('Time t')

% A solution profile can also be illuminating.
figure
plot(x,u(end,:))
title('Solution at t = 2')
xlabel('Distance x')
ylabel('u(x,2)')
% --------------------------------------------------------------
function [c,f,s] = pdex1pde(x,t,u,DuDx)
v=0.2;
lambda=0.5;
c = 1;
f = -v*u;
s = -lambda*u+lambda/v*u;
% --------------------------------------------------------------
function u0 = pdex1ic(x)
u0 = (sin(pi*2*x))*(x<0.5)*(x>0); %(x<0.51)*(x>0.49);
% --------------------------------------------------------------
function [pl,ql,pr,qr] = pdex1bc(xl,ul,xr,ur,t)
pl = 0;
ql = 1;
pr = ur;
qr = 0;
