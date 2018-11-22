Summary of the key parts of this directory relevant for thesis chapter on fine grained model.

Parts of this directory can be used/have been used to implement a version of lazy abc for the 3D version of the fine grained VJP model.
This can be used to infer parameters for this model from synthetic data, and show the effects of different summary statistics.

%%%%%%%%%%%%%%%%%%%%%%
Notes
Unclear: inference is for 3 parameters. The full model seems to rely on ~7? parameters. Which parameters are fixed and how is this done?
Unclear: priors are uniform on log of parameters. Would more informative priors circumvent the need for using lazy ABC?

For thesis: want to run lazy ABC and compare different summary statistics on synthetic and real data.
Can achieve this by using combined summary statistic results and then doing some post processing.

Generate data (Synthetic or real)
Sample from prior. (3 params. 7 params?)
Run model. (3 params, 7 params?)
Calculate summary statistics.
Store these.
Restore from save.
Select closest subset (alpha=?,summary statistics) 

%%%%%%%%%%%%%%%%%%%%%%

lazy_abc.m - main script to run ABC in batches in parallel, calls other scripts
propose_theta_lazy.m - proposes new parameters, calculates typical step size which is used to determine whether to continue a simulation
continue_lazy.m - if you have continued to the simulation, this will call the model to create the data and get the distance from observed 
create_synthetic_data.m - wrapper that calls the VJP model and calculates summary statistics
my_dist.m - calculates euclidean distance between data sets
../summary_statistic_calculate_combined_3Dv2.m - call the VJP model to calculate summary statistics
posterior_contour_plot.m - for ploting the posterior distribution
 
