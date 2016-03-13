function [dist] = continue_lazy(params,data_summary,opts)

%generate data
q = create_synthetic_data(params,opts);
dist = my_dist(q,data_summary);
end
