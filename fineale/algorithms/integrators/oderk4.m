% Fixed-step fourth-order Runge-Kutta  integrator.
%
% Function to integrate a system of first-order
% initial-value ODE's by the fourth-order Runge-Kutta  method.
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
function [ts, ys] = oderk4 (rhsfun, tspan, y0, options, varargin)
    f=rhsfun;
    try,        h=options.(matchfieldname(options, 'InitialStep'));
    catch,      error(['Option InitialStep must be supplied']);    end
    % check whether we are going forward or backward in time
    if (tspan(2)<tspan(1)),h=-abs(h);end
    nsteps = abs(ceil ((tspan(2) - tspan(1)) / h));
    ts = zeros(nsteps+1,1);
    ys = zeros(nsteps+1,length(y0));
    t = tspan(1);
    ts(1) = t;
    ys(1,:) = y0';
    for step = 2:nsteps+1
        if (abs(tspan(2)-ts(step-1))<abs(h)),h=tspan(2)-ts(step-1);end
        ts(step) = ts(step-1) + h;
        ys(step,:) = rk4step(f, h, ts(step-1), ys(step-1,:)', varargin{:});
    end
end

function [yn1] = rk4step (rhsfun, h, tn, yn, varargin)
    k1 = feval(rhsfun, tn, yn);
    k2 = feval(rhsfun, tn+h/2, yn+h*k1/2, varargin{:});
    k3 = feval(rhsfun, tn+h/2, yn+h*k2/2, varargin{:});
    k4 = feval(rhsfun, tn+h, yn+h*k3, varargin{:});
    yn1 = yn + h/6*(k1+2*k2+2*k3+k4);
end