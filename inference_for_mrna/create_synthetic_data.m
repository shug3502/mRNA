function data_summary = create_synthetic_data(params,opts)

%defaults for calculating summary stat:
% for par_params expect a vector of params of length 7
%is_parallel=1;
%option_a_summary_statistic=1; %use distribution at various times as summary stat
    %generate summary stat
    q = summary_statistic_calculator_combined_3Dv2(params,opts.num_particles,opts.is_parallel);
%   q = summary_statistic_calculator_3D(params,opts.num_particles,opts.is_parallel,opts.ss);
    data_summary = reshape(q,1,[]);
end

