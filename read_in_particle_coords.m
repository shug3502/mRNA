function [q,x,y] = read_in_particle_coords(file_num,pixel_sizes)
%created 10/6/2015
%last edit 11/6/2015
% reads in centroids data ouput by python and uses this as positions of
% particles

%Issues with this: where to define axis? Where is the single cell that we
%are interested in? 
%Also which time point does the given distribution correspond to and how
%can we obtain this?

%Requires input of pixel sizes


x_scale = pixel_sizes(1); %x scale is how many pixels correspond to a micron?
y_scale = pixel_sizes(2);

fname = sprintf('Data/my_mRNA_output%d.csv',file_num);
M = csvread(fname,1,0);
fclose('all');

y = M(:,2)/x_scale;
x = M(:,3)/y_scale;

%split into bins
%params for splitting into bins
delx = 1;
L= 52;
if length(x)/pixel_sizes(1)<=0 || length(x)/pixel_sizes(1)>L
    error('min x <0 or max x > L');
end

[Num_in_bins,~] = histc(x,0:delx:L);
q = Num_in_bins/sum(Num_in_bins)/delx; %estimate of q at time T
% figure
% bar(q)
end