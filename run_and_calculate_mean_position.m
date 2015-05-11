function run_and_calulate_mean_position(N,num_modes)
tic;
close all
dt = 0.02;
time_vec = 0:dt:0.5;
%parameters
params.nu1 = 0.4;
params.nu2 = 0.08;
params.lambda = 1/0.13;
params.omega = 1/6;
params.x_0=0.5;
params.theta_0 = 0;
params.phi = 0.58 ; %bias in MTs

L=30; %length of cell
delx = 1; %bin size for integral

mean_position = zeros(size(time_vec));
second_moment = zeros(size(time_vec));
q_estimate = zeros(L+1,length(time_vec));

for k = 1:length(time_vec)
    final_xpos = zeros(N,1);
    final_ypos = zeros(N,1);
    for j=1:N
        [~, ~, final_xpos(j), final_ypos(j), pathx, jump_times] = velocityjump2D_ModesInput(time_vec(k), params, 1, 0, num_modes, 0, 3*j);
    end
    
%     figure(1);
%     subplot(length(time_vec),1,k)
%     hist(final_xpos,0:5*delx:L);
    
    %split into bins
    [Num_in_bins,edges] = histc(final_xpos,0:delx:L);
    q_estimate(:,k) = Num_in_bins/N; %estimate of q at time T
    mean_position(k) = (delx/2:delx:(L-delx/2))*q_estimate(1:L,k)*delx+delx*L*q_estimate(L+1,k); %use centre of each bin 
    second_moment(k) = (delx/2:delx:(L-delx/2)).^2*q_estimate(1:L,k)*delx+delx*L^2*q_estimate(L+1,k); %use centre of each bin    
end

speed_estimate = diff(mean_position)/dt
%sum(speed_estimate(1:8))/8
figure(2)
t = time_vec*60^2;
subplot(2,1,1);
errorbar(t,mean_position,sqrt(second_moment-mean_position.^2),'linewidth',3)
%plot(time_vec,mean_position,'r','linewidth',3)
grid on
hold on

%PLOT ANALYTIC RESULT
mu_initial = params.x_0; %initial mean position
nu_1 = params.nu1; %speed in active transport mode
F1 = (4*params.phi-2)/pi;
mu1 = min((nu_1*F1*t + mu_initial),L); %note time scaled to seconds
if num_modes>1
    mu1 = min((0.5*nu_1*F1*t + mu_initial),L); %adjust for different analytical results for different numbers of modes
    mud = min((nu_1*F1*t + mu_initial),L);
    plot(t,mud,'m--','linewidth',3);
end
hold on
plot(t, mu1, 'g--', 'linewidth',3);
set(gca, 'fontsize',14)
xlabel('t');
ylabel('mean position');

for s=1:length(time_vec)
subplot(2,1,2);
bar(0:delx:L,q_estimate(:,s))
set(gca, 'fontsize',14)
xlabel('x')
ylabel('q')
grid on
axis([0,35,0,1])
pause(0.4)
end

figure(3)
plot(t,second_moment,'r','linewidth',3)
grid on
hold on

%PLOT ANALYTIC RESULT
mu_initial = params.x_0; %initial mean position
nu_1 = params.nu1; %speed in active transport mode
omega = params.omega; %switching rate
F1 = (4*params.phi-2)/pi;
mu2 = min((nu_1*F1*t).^2 + 2*nu_1*F1*t*mu_initial+ mu_initial^2,L^2); %note time scaled to seconds
if num_modes>1
   %adjust for different analytical results for different numbers of modes

    mu2_prime = min(( mu_initial^2 + nu_1*F1*( mu_initial*t ... 
     + nu_1/2*F1*(t).^2 ...
      + mu_initial/(2*omega)*(1-exp(-2*omega*t)) ... 
      + t/(2*omega) + 1/(4*omega^2)*(exp(-2*omega*t)))), L^2); 
%      + nu_1*(4*phi/100-2)/pi*time_vec*60^2.*exp(2*omega*time_vec*60^2)/(2*omega) ...
%      - mu_initial/(2*omega) - nu_1*(4*phi/100-2)/pi*exp(2*omega*time_vec*60^2)/(4*omega^2)...
%      + nu_1/(4*omega^2)*(4*phi/100-2)/pi),L^2);
plot(t, mu2_prime, 'k--', 'linewidth',3);

mu2 = min(mu_initial^2 + nu_1/2*F1*( 2*mu_initial*t + nu_1/2*F1*t.^2),L^2);

end

plot(t, mu2, 'g-.', 'linewidth',3);
set(gca, 'fontsize',14)
xlabel('t');
ylabel('second moment');

second_moment
mu2
toc;
end


