function dgdt = ode_coarse_mRNA(t,g,params)
%% Created 11/8/16 JH
%% last edit 11/8/16
%% coarse grained model using odes
%% assume fast transport?
%% what is relative contribution to rate of accumulation from transport and from production?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<1
    t = 0;
    g = ones(16,1);
    params.a = 1;
    params.b = 1;
    params.d = 0;
end

A = ones(16,1); A(1)=0; %matrix of which cells produce RNA
B = get_nc_transitions(0)'; %assume cant exit oocyte
%B(1,:) = zeros(1,16);
D = ones(16,1);
% dgdt = params.a*A + params.b*B*g - params.d*D.*g; %dimensional equation

dgdt = params.a/params.b*A + B*g - params.d*D.*g; %nondimensional equation



