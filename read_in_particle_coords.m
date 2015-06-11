function [x,y] = read_in_particle_coords(pixel_sizes)

x_scale = pixel_sizes(1);
y_scale = pixel_sizes(2);

    fname = sprintf('my_mRNA_output1.csv');
    M = csvread(fname,1,0);
    fclose('all');
    
    y = M(:,2)/x_scale;
    x = M(:,3)/y_scale;


end