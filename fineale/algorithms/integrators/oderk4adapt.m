% Adaptive fourth-order Runge-Kutta  integrator.
%
% Function to integrate a system of first-order
% initial-value ODE's by the fourth-order Runge-Kutta  method.
%
% The integrator accepts the same arguments as the built-in
% Matlab solvers ode23, ...  The time step needs to be set
% by the options argument, field InitialStep. The absolute tolerance needs
% to be set by the options argument, field AbsTol. Only the norm of the
% solution is controlled by AbsTol.
% 
% This is a time-step-adaptive integrator: the time step is adjusted
% according to the formula
%      dt_estim = dt*((1-1/2^4)*AbsTol/norm(tyn1-yn1))^(1/5);
% where 
% yn1 = solution predicted with the full step dt,
% tyn1 = solution predicted with the two half steps dt/2,
% 
% The solution accepted for the next step is tyn1.
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
function [ts, ys] = oderk4adapt (rhsfun, tspan, y0, options, varargin)
    f=rhsfun;
    try,        h=options.(matchfieldname(options, 'InitialStep'));
    catch,      error(['Option InitialStep must be supplied']);    end
    hmax = inf;
    if (isfield(options,'MaxStep'))
        hmax= options.MaxStep;
    end
    AbsTol = 1e-3;
    try,        AbsTol = options.(matchfieldname(options, 'AbsTol'));
    catch,    AbsTol = 1e-3;     end
    if (AbsTol<=0), AbsTol = 1e-3;     end
    nsteps = ceil ((tspan(2) - tspan(1)) / h);
    ts = zeros(nsteps+1,1);
    ys = zeros(nsteps+1,length(y0));
    t = tspan(1);
    ts(1) = t;
    ys(1,:) = y0';
    step =2;
    while true
        tn=ts(step-1);
        yn =ys(step-1,:)';
        yn1 =rk4step (f, h, tn, yn, varargin{:});
        tyn12 =rk4step (f, h/2, tn, yn, varargin{:});
        tyn1 =rk4step (f, h/2, tn+h/2, tyn12, varargin{:});
        hh =h*((1-1/2^4)*AbsTol/norm(tyn1-yn1))^(1/5);
        if hh>=h % accept this step
            ts(step) = tn+ h;
            ys(step,:) = tyn1';
            step=step+1;
        end
        h=min([hh,2*h]);
        if (h>hmax)
            h = hmax;
        end
        if (ts(step-1)+h > tspan(2))
            h=tspan(2)-ts(step-1);
        end
        if (ts(step-1) == tspan(2))
            ts=ts(1:step-1);
            ys =ys(1:step-1,:);
            break;
        end
    end
end
%
% Do one 4th-order Runge-Kutta step for the equation dy/dt = f(t,y).
% 
% y_n+1 = y_n + h/6*(k1+2*k2+2*k3+k4);
%
% Arguments: 
%   rhsfun = rhs function f(t,y)
%   h      = step length
%   tn     = time t_n
%   yn     = function value y_n
%   
function [yn1] = rk4step (rhsfun, h, tn, yn, varargin)
    k1 = feval(rhsfun, tn, yn);
    k2 = feval(rhsfun, tn+h/2, yn+h*k1/2, varargin{:});
    k3 = feval(rhsfun, tn+h/2, yn+h*k2/2, varargin{:});
    k4 = feval(rhsfun, tn+h, yn+h*k3, varargin{:});
    yn1 = yn + h/6*(k1+2*k2+2*k3+k4);
end
