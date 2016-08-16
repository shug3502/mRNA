function perfectly_processed_pdfs(fig_name)
%% From show and tell coding...
%% Chris's Decimal Delight
%
set(gca,'XTickLabel',arrayfun(@(s)sprintf('%0.1f', s), cellfun(@(s)str2num(s), get(gca,'XTickLabel')), 'UniformOutput', false))

%% Ruth's Perfectly Processed PDFs
% How to output pdfs from Matlab that are the same size etc
%set(gca,'LooseInset',get(gca,'TightInset'))
%set(gca,'FontSize',20)
%set(gcf,'PaperUnits','centimeters','PaperSize',[20/3+0.2 16/3+0.2], 'PaperPosition', [0.1 0.1 20/3+0.1 16/3+0.1]);
saveas(gcf,sprintf('%s.pdf',fig_name),'pdf');
%%
print(sprintf('%s.eps',fig_name),'-depsc');
