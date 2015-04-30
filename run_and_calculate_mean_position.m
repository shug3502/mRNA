function run_and_calulate_mean_position(N)
tic;
close all
dt = 0.02;
time_vec = 0:dt:0.3;
%parameters
phi = 58; %bias in MTs
L=30; %length of cell
delx = 1; %bin size for integral

mean_position = zeros(size(time_vec));
q_estimate = zeros(L+1,length(time_vec));

for k = 1:length(time_vec)
    final_xpos = zeros(N,1);
    final_ypos = zeros(N,1);
    for j=1:N
        [~, ~, final_xpos(j), final_ypos(j)] = velocityjump2D_ModesInput(time_vec(k), phi, 1, 0, 2, 0, 3*j);
    end
    
%     figure(1);
%     subplot(length(time_vec),1,k)
%     hist(final_xpos,0:5*delx:L);
    
    %split into bins
    [Num_in_bins,edges] = histc(final_xpos,0:delx:L);
    q_estimate(:,k) = Num_in_bins/N; %estimate of q at time T
    mean_position(k) = (delx/2:delx:(L-delx/2))*q_estimate(1:L,k)*delx+delx*L*q_estimate(L+1,k); %use centre of each bin 
end
speed_estimate = diff(mean_position)/dt
sum(speed_estimate(1:8))/8
figure(2)
subplot(2,1,1);
plot(time_vec,mean_position,'r','linewidth',3)
xlabel('t')
ylabel('mean position')
grid on
for s=1:length(time_vec)
subplot(2,1,2);
bar(0:delx:L,q_estimate(:,s))
xlabel('x')
ylabel('q')
grid on
axis([0,35,0,1])
pause(0.4)
end
toc;
end


