function solve_mfpt_pde
%created JH 30/4/15
%last edit 30/4/15
%solving evolution equation avergaed over theta
%ie just transport equation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L=30; %30 microns
dy = 0.1;
num_alph_pts = 101
alphas = linspace(0,2*pi,num_alph_pts);
num_space_pts = L/dy+1
lambda = 1/0.13;

phi = 58;
speed = 0.4; %microns per second
F1 = (4*(phi/100)-2)/pi; %integral on sin(theta)*T(theta) d theta 
F2 = pi*(3/2 - phi); %integral of theta *T(theta) d theta

tau = zeros(num_space_pts,num_alph_pts);

%tau(:,1) = poisspdf(1:num_space_pts,2/dy);
tau(1,:) = 100*sin(alphas);

for j=2:num_space_pts
% size(G2(alphas,lambda,F2))
% size(G1(alphas,speed,lambda,F1))
% size(tau(j-1,1:(num_alph_pts-2)))
tau(j,2:(num_alph_pts-1)) = 0.5*(tau(j-1,1:(num_alph_pts-2))+tau(j-1,1:(num_alph_pts-2))) ...
    - dy./G1(alphas(2:end-1),speed,lambda,F1).*(1+G2(alphas(2:end-1),lambda,F2).*...
    (0.5*(tau(j-1,1:(num_alph_pts-2)) - tau(j-1,3:(num_alph_pts)))));
tau(num_space_pts,:) = 0;
tau(:,1) = tau(:,num_alph_pts);

end
close all
figure;
plot(0:dy:L,tau(:,1),'g-') %,0:dy:L,tau(:,num_alph_pts))
xlabel('x')
ylabel('probability')

end

function g = G1(alpha,speed,lambda,F1)
g = speed*(sin(alpha)*(1-lambda)+lambda*(F1));
end

function g = G2(alpha,lambda,F2)
g = lambda*(F2 - alpha);
end
