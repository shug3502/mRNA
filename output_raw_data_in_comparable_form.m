function output_raw_data_in_comparable_form

[q,x,y] = read_in_particle_coords(9,[10,10]);

figure
imagesc(1,0:52,log(q)) % plot whole time course
colormap(gray)
set(gca, 'fontsize', 20);
set(gca,'XTick',25)
yticklabels = 0:10:52;
yticks = linspace(1, size(q, 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)))
xlabel('t=25')
ylabel('Position')

figure
imagesc(1,0:52,(q)) % plot whole time course
colormap(gray)
set(gca, 'fontsize', 20);
set(gca,'XTick',25)
yticklabels = 0:10:52;
yticks = linspace(1, size(q, 1), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', flipud(yticklabels(:)))
 xlabel('t=25')
 ylabel('Position')


