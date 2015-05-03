%attempt 2 at pdepe
function attempt2_at_pdepe
m=0;
dx=0.5;
dt = 0.1;
L=30;
T=10;
xmesh = 0:dx:L;
tspan = 0:dt:T;

sol=pdepe(m,@pdefun,@icfun,@bcfun,xmesh,tspan);

sol
surf(xmesh,tspan,sol)

function [mu,f,s]=pdefun(x,t,u,ux)
mu = 1;
speed = 0.4;
F1 = 0.01;
f = -speed*F1*u;
s = 0;
end

function phi=icfun(x)
phi = sin(x)*(x>5)*(x<10);
end

function [pa,qa,pb,qb]=bcfun(a,ua,b,ub,t)
pa = 0;
qa = ua;
pb = 1;
qb = 0;
end
end