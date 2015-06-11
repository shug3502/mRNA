function [q,x,y] = read_in_particle_coords(pixel_sizes)

x_scale = pixel_sizes(1); %x scale is how many pixels correspond to a micron?
y_scale = pixel_sizes(2);

fname = sprintf('my_mRNA_output1.csv');
M = csvread(fname,1,0);
fclose('all');

y = M(:,2)/x_scale;
x = M(:,3)/y_scale;

%split into bins
%params for splitting into bins
delx = 1;
L= 52;
if length(x)/pixel_sizes(1)<=0 | length(x)/pixel_sizes(1)>L
    error('min x <0 or max x > L');
end

[Num_in_bins,~] = histc(x,0:delx:L);
q = Num_in_bins/sum(Num_in_bins)/delx; %estimate of q at time T
figure
bar(q)
end