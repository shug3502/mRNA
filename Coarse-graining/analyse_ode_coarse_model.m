%% analyse_ode_coarse_model
%% Created 11/8/16 JH
%% last edit 11/8/16
%% analyse and integrate coarse grained model using odes

%variables to change: initial condition, final time of simulation, a/b

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g0 = ones(16,1); %initially no/some grk mRNA present
params.a = 0.01;
params.b = 1;
params.d = 0;
opts=[]; %odeset('NonNegative',1);
[t,g]=ode45(@ode_coarse_mRNA,[0,10],g0,opts,params);

%find eigenvalues for general solution
A = ones(16,1); A(1)=0;
B = get_nc_transitions(0)'; %assume cant exit oocyte
%B(1,:) = zeros(1,16);

[V,D] = eig(B);
k = -B\(params.a/params.b*A); %note B is singular
k(1) = 0; %since B is singular (k(1) gets multiplied by 0), this is arbitrary
C = V\(g0-k);

%plot rate of increase of rna in oocyte
dgdt=zeros(size(g));
my_soln = zeros(size(g));
for j=1:size(g,1)
dgdt(j,:) = ode_coarse_mRNA(t(j),g(j,:)',params);
my_soln(j,:) = (V*diag(exp(diag(D*t(j)))))*C + k + [sum(A)*params.a/params.b; zeros(15,1)]*t(j);
end

my_soln-g;

close all;
figure;
plot(t,g(:,1),'linewidth',3);
hold all
plot(t,my_soln(:,1),'linewidth',3);
%plot(t,k(2)*ones(size(t)),'k--','linewidth',3);
plot(t,dgdt(:,1),'linewidth',3);
box on
set(gca,'fontsize',20);
xlabel('t');
ylabel('y');
legend('g numeric','g analytic','dg/dt','Location','NorthWest');
title(sprintf('a/b = %f', params.a/params.b));

