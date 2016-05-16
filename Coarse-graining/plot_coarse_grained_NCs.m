function plot_coarse_grained_NCs(save_name, omit_oocyte)
%created 11/5/16
%last edited 11/5/16
%JH
%plots data for distn of mRNA in nurse cells onto a representation of NCs
%see caceres development paper 2005 for NCs map
%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;

load(save_name); %load previous data including params etc

%map integers 1:16 onto their position in NCs layout diagram
NCs_map = [2,3;
    2,2;
    3,3;
    3,2;
    2,4;
    2,1;
    3,4;
    3,1;
    1,3;
    1,2;
    4,3;
    4,2;
    1,4;
    1,1;
    4,4;
    4,1];

M = zeros(4); %store mRNA concentration (proportion of total mRNA) in M
for i=1:size(distn,1)
    %loop over the cells
    for j=1:16
        M(NCs_map(j,1),NCs_map(j,2)) = distn(i,j);
    end
    %figure(i+size(distn,1));
    %imagesc(M);
    % or...
    figure(i);
    set(gca,'XTickLabel','');
    set(gca,'YTickLabel','');
    title(sprintf('Nurse cells at t=%2.1f ',params.dt*params.thin_results*(i-1)));
    hold on;
    if omit_oocyte;
        M(NCs_map(1,1),NCs_map(1,2))=0; %remove oocyte for plotting
        mm=max(max(M));
    else
        mm=max(max(distn));
    end
    for i1=1:4
        for i2=1:4
            rectangle('Position',[i1-1,i2-1,1,1],...
                'Curvature',[1,1], 'FaceColor',M(5-i2,i1)/mm*ones(1,3)); %plots/counts from bottom left, rather than top left
        end
    end
    
    %     figure;
    %     rectangle('Position',[NCs_map(j,1)-1,NCs_map(j,2)-1,1,1],...
    %   'Curvature',[1,1], 'FaceColor','r')
end