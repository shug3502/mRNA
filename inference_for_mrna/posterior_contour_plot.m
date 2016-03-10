function posterior_contour_plot(posterior, save_name)
%created 11/12/16
%last edit 11/2/16

%take posterior samples eg. from ABC and plot kernel estimate
%should be adjusted to account for particle weights

%%%%%%%%%%%%%%%%%%%%%

figure;    
indices = [1,2; 1,3; 2,3];
labels = cell(3);
labels{1} = '\nu';
labels{2} = '\lambda';
labels{3} = '\phi';

for i=1:3
    data= posterior(:,indices(i,:));
    % call the routine
    [bandwidth,density,X,Y]=kde2d(data);
    % plot the data and the density estimate
    subplot(3,1,i);
    contour3(X,Y,density,50), hold on
    %plot(data(:,1),data(:,2),'r.','MarkerSize',5)    
    xlabel(labels{indices(i,1)});
    ylabel(labels{indices(i,2)});
    view(2);
    colorbar;
end
if nargin>1
	print(save_name,'-depsc');
end
