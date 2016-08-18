%% inverse problem to find a and b
%% Created 16/8/16 JH
%% Last edit 16/8/16
%% want to find a and b to get ratio a/b

%%approach empirically: find how a summary statistic depends on b
%%%%%%%%%%%%%%%%%%

n=1;
b=10.^linspace(0,0,n);
T_half = zeros(n,1);
for j=1:n %loop over values of b
%generate synthetic data
g0 = zeros(16,1); g0(4)=1; %injection expmt
real.a = 0;
real.b = b(j);
real.d = 0;
opts=[]; %odeset('NonNegative',1);
[t,g]=ode45(@ode_coarse_mRNA,[0,100],g0,opts,real);
%corrupt with some noise
noise_sd = 0.0;
g = g.*exp(noise_sd*randn(size(g)));

ind_half = find(g(:,1)>0.5,1);
T_half(j) = t(ind_half)/real.b;
end
%visualise
close all;
figure;
plot(t/real.b,g(:,1),'linewidth',3);
hold all
box on
set(gca,'fontsize',20);
xlabel('t');
ylabel('y');

figure;
plot(log10(b),log10(T_half),'bo','markersize',12);
hold all
box on
set(gca,'fontsize',20);
xlabel('$\log_{10} b$','interpreter','LaTex');
ylabel('$\log_{10} T_{\frac{1}{2}}$','interpreter','LaTex');

