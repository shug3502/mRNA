function numerically_int_q_guess(N,num_modes)

%have made a sigmoidal ansatz for q(x,t) = 1/(1+exp(0.5*(x-nu*F1*t)))
L=30;
delx = 0.01;
params.nu1 = 0.4;
params.nu2 = 0.08;
params.lambda = 1/0.13;
params.omega = 1/6;
params.x_0=0.0005;
params.theta_0 = 0;
params.phi = 0.58 ; %bias in MTs

F1 = (4*params.phi-2)/pi;

beta = 0.5;
x = 0:delx:L;
t = (0:0.02:0.5)*60^2;


q = 1./(1+exp(beta*(repmat(x,length(t),1)-params.nu1*F1*repmat(t',1,length(x)))));

close all
figure;
plot(x,q(10,:),'linewidth',3)
grid on

figure;
plot(t,sum(repmat(x,length(t),1).*q*delx,2)./sum(q*delx,2),'g--','linewidth',3)
grid on
hold on

dt = 0.02;
time_vec = 0:dt:0.5;
t = time_vec*60^2;
l_t = length(time_vec);
%parameters
params.nu1 = 0.4;
params.nu2 = 0.08;
params.lambda = 1/0.13;
params.omega = 1/6;
params.x_0=0.0005;
params.theta_0 = 0;
params.phi = 0.58 ; %bias in MTs

%L=30; %length of cell
delx = 0.1; %bin size for integral

mean_position = zeros(1,l_t);
second_moment = zeros(1,l_t);
jumps = zeros(N,10^4);
q_estimate = zeros(L/delx+1,l_t);
xpos = zeros(N,10^4);
xpos_discrete_time = zeros(N,10^4);

for j=1:N
    [~, ~, ~, ~, xpos_temp, jump_temp] = velocityjump2D_ModesInput(max(time_vec), params, 1, 0, num_modes, 0, 13*j);
    jumps(j,1:length(jump_temp)) = jump_temp;
    xpos(j,1:length(xpos_temp)) = xpos_temp;    
    for w=2:l_t
        if isempty(xpos(j,find(t(w)<jumps(j,:),1,'first')))
            xpos_discrete_time(j,w) = L;
        else
            xpos_discrete_time(j,w) =  xpos(j,find(t(w)<jumps(j,:),1,'first'));
        end
    end
end

for w=1:l_t
    %split into bins
    [Num_in_bins,edges] = histc(xpos_discrete_time(:,w),0:delx:L);
    q_estimate(:,w) = Num_in_bins/N/delx; %estimate of q at time T

    mean_position(w) = (delx/2:delx:(L-delx/2))*q_estimate(1:(end-1),w)*delx+delx*L*q_estimate(end,w); %use centre of each bin
    second_moment(w) = (delx/2:delx:(L-delx/2)).^2*q_estimate(1:(end-1),w)*delx+delx*L^2*q_estimate(end,w); %use centre of each bin
end

t = time_vec*60^2;
errorbar(t,mean_position,sqrt(second_moment-mean_position.^2),'linewidth',3)

end
