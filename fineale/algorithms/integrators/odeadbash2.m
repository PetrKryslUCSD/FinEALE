% Second-order Adams-Bashforth method.
%
% Function to integrate a system of first-order 
% initial-value ODE's by the second-order Adams-Bashforth method.
% 
% One second-order Adams-Bashforth step for the equation dy/dt = f(t,y) is
% 
% y_n+1 = y_n + 3/2*h * f(t_n, y_n) - 1/2*h f(t_n-1, y_n-1)
%
% The first step is done from slope calculated with one forward Euler step
% backwards.
% 
% The integrator accepts the same arguments as the built-in
% Matlab solvers ode23,….  The time step needs to be set
% by the options argument, field InitialStep.
% 
% Arguments:
%    rhsfun  = right-hand side function: returns f(t,y) 
%              as a column vector.
%    tspan   = time span (start time, usually zero, 
%              and stop time).
%    y0      = initial condition.
%    h       = step size; constant.
%
%
% Copyright (C) 2007, Petr Krysl
% 
function [ts, ys] = odeadbash2 (rhsfun, tspan, y0, options, varargin)
    h = options.InitialStep(1);
    nsteps = ceil ((tspan(2) - tspan(1)) / h);
    ts = zeros(nsteps+1,1);
    ys = zeros(nsteps+1,length(y0));
    t = tspan(1);
    step = 1;
    ts(step) = t;
    ys(step,:) = y0';
    % start by estimating the slope at the time t_n-1 by a forward Euler
    % step (well, backward)
    fn = feval (rhsfun, ts(step), ys(step,:)', varargin{:});
    fnm1 = feval (rhsfun, ts(step)-h, ys(step,:)'-h*fn, varargin{:});
    for step = 2:nsteps+1
        ts(step) = ts(step-1) + h;
        ys(step,:) = (feval ('AdBash2', rhsfun, h, ts(step-1), ys(step-1,:)', fn, fnm1, varargin{:}))';
        fnm1 = fn;
        fn = feval (rhsfun, ts(step), ys(step,:)', varargin{:});
    end
    return;
end

%
% Do one  second-order Adams-Bashforth step for the equation dy/dt = f(t,y).
% 
% y_n+1 = y_n + 3/2*h * f(t_n, y_n) - 1/2*h f(t_n-1, y_n-1)
%
% Arguments: 
%   rhsfun = rhs function f(t,y)
%   h      = step length
%   tn     = time t_n
%   yn     = function value y_n
%   
function [yn1] = AdBash2 (rhsfun, h, tn, yn, fn, fnm1, varargin)
    yn1 = yn + (3/2)*h * fn - (1/2)*h * fnm1;
end
