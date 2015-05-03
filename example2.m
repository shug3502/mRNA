function example2 ( )
timestamp();
  m = 0;
%
%  Define the spatial mesh.
%
  nx = 1001;
  L = 30; % microns;
  xmesh = linspace ( 0, L, nx );
%
%  Define the time mesh.
%
  nt = 101;
  tspan = linspace ( 0.0, 1, nt );
%
%  Call PDEPE() for the solution.
%
  sol = pdepe ( m, @pdefun, @icfun, @bcfun, xmesh, tspan );
%
%  Even though our solution is "really" a 2D array, PDEPE stores it
%  in a 3D array SOL(:,:,:).  The surf() command needs a 2D array to plot,
%  so let's copy U out of SOL.
%
  u = sol(:,:,1);

  figure ( 1 )
  surf ( xmesh, tspan, u, 'EdgeColor', 'None' );
  title ( 'Example 2: Solution Over Time', 'Fontsize', 16 );
  xlabel ( '<--- X --->' )
  ylabel ( '<--- T --->' );
  zlabel ( '<---U(X,T)--->' );
  filename = 'example2.png';
  print ( '-dpng', filename );
  fprintf ( 1, '\n' );
  fprintf ( 1, '  Saved solution plot in file "%s"\n', filename );
%
%  Plot the initial condition, U at time 0.
%
  figure ( 2 )
  plot ( xmesh, u(1,:), 'LineWidth', 3 );
  grid on
  title ( 'Example 2: Initial Condition', 'Fontsize', 16 );
  xlabel ( '<--- X --->' )
  ylabel ( '<--- U(X,T0) --->' );
  filename = 'example2_ic.png';
  print ( '-dpng', filename );
  fprintf ( 1, '  Saved initial condition plot in file "%s"\n', filename );
%
%  Plot the solution U at a fixed point, with time varying.
%
  
%   figure ( 3 )
%   mid = 1 + floor ( 55 * ( nx - 1 ) / 100 );
%   plot ( tspan, u(:,mid), 'LineWidth', 3 );
%   grid on
%   title ( 'Example 2: Time evolution of solution at X=5.0', 'Fontsize', 16 );
%   xlabel ( '<--- T --->' )
%   ylabel ( '<--- U(5.0,T) --->' );
%   filename = 'example2_profile.png';
%   print ( '-dpng', filename );
%   fprintf ( 1, '  Saved time evolution plot in file "%s"\n', filename );
%
%  Animate the profile.
%  I wish I could also display the running value of time, but it
%  does not seem possible.
%
  figure ( 4 )
  fig = plot ( xmesh, u(1,:), 'erase', 'xor' );
  title ( 'Profile animation', 'Fontsize', 16 );
  grid on
  for k = 2 : length ( tspan )
    set ( fig, 'xdata', xmesh, 'ydata', u(k,:) );
    pause ( 0.1 );
  end
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'EXAMPLE2:\n' );%
  fprintf ( 1, '  Normal end of execution.\n' );
  fprintf ( 1, '\n' );
  timestamp ( );

  return
end
function [ c, f, s ] = pdefun ( x, t, u, dudx )

%*****************************************************************************80
%
%% PDEFUN defines the components of the PDE.
%
%  Discussion:
%
%    The PDE has the form:
%
%      c * du/dt = x^(-m) d/dx ( x^m f ) + s
%
%    where m is 0, 1 or 2,
%    c, f and s are functions of x, t, u, and dudx, 
%    and most typically, f = dudx.
%
  c = 1.0;
  f = - 0.04 * u;
  s = 0.0;

  return
end
function u0 = icfun ( x )

%*****************************************************************************80
%
%% ICFUN defines the initial conditions.

  %u0 = 1.0 / ( 1.0 + ( x - 5.0 ).^2 );
u0 = (1.0 / ( 1.0 + ( x - 5.0 ).^2 )).*(x>5.0)*(x<10.0);
  
  return
end
function [ pl, ql, pr, qr ] = bcfun ( xl, ul, xr, ur, t )

%*****************************************************************************80
%
%% BCFUN defines the boundary conditions.

  pl = 0.0; %ul
  ql = ul; %0.0
  pr = ur;
  qr = 0.0;

  return
end
function timestamp ( )

%*****************************************************************************80
%
%% TIMESTAMP prints the current YMDHMS date as a timestamp.
  t = now;
  c = datevec ( t );
  s = datestr ( c, 0 );
  fprintf ( 1, '%s\n', s );

  return
end
