function D = my_dist(data1,data2)
scale=1;
D = norm((data1-data2)./scale,2);
end

