% Modified Euler integrator.   Alternative formulation.
%
% Function to integrate a system of first-order 
% initial-value ODE's by the modified Euler method (yet another one).
% 
% One modified Euler step for the equation dy/dt = f(t,y) is
% 
%    y_n+1 = y_n + h * f(t_n+1/2, y_n+1/2) 
% 
% which is approximated as 
%
%    y_n+1 = y_n + h * f(t_n+1/2, y_p) 
% 
% where y_p = y_n + h/2 * f(t_n, y_n)
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
% Copyright (C) 2009, Petr Krysl
% 
function [ts, ys] = odemeul2 (rhsfun, tspan, y0, options, varargin)
    try,        h=options.(matchfieldname(options, 'InitialStep'));
    catch,      error(['Option InitialStep must be supplied']);    end
    nsteps = ceil ((tspan(2) - tspan(1)) / h);
    ts = zeros(nsteps+1,1);
    ys = zeros(nsteps+1,length(y0));
    t = tspan(1);
    ts(1) = t;
    ys(1,:) = y0';
    for step = 2:nsteps+1
        if (ts(step-1) + h >tspan(2)), h=tspan(2)-ts(step-1); end
        ts(step) = ts(step-1) + h;
        ys(step,:) = (MEuler2(rhsfun, h, ts(step-1), ys(step-1,:)', varargin{:}))';
    end
    return;
end

%
% Do one Modified Euler step for the equation dy/dt = f(t,y). (Yet another
% method that got the name of the modified Euler...)
% 
% y_n+1 = y_n + h * f(t_n+1/2, y_n+1/2) 
% 
% which is approximated as 
%
% y_n+1 = y_n + h * f(t_n+1/2, y_p) 
% 
% where y_p = y_n + h/2 * f(t_n, y_n)
%
% Arguments: 
%   rhsfun = rhs function f(t,y)
%   h      = step length
%   tn     = time t_n
%   yn     = function value y_n
%   
function [yn1] = MEuler2 (rhsfun, h, tn, yn, varargin)
    tp  = tn + h/2;
    fn  = feval (rhsfun, tn, yn, varargin{:});
    yp  = yn + h/2 * fn;
    yn1 = yn + h * feval (rhsfun, tp, yp, varargin{:});
end
