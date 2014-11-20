% Trapezoidal-rule time integrator.
%
% Function to integrate a system of first-order 
% initial-value ODE's by the trapezoidal  method.
% 
% One trapezoidal method step for the equation dy/dt = f(t,y) is
% 
% y_n+1 = y_n + h/2 * (f(t_n+1, y_n+1) +f(t_n, y_n))
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
%    options = as returned by odeset(). Field InitialStep needs to be set
%              to the desired step size which kept throughout the
%              calculation.  The field  RelTol may be set to control the
%              relative tolerance of the solution components (1e-3 is the
%              default). 
%
% Copyright (C) 2007, Petr Krysl
% 
function [ts, ys] = odetrap(rhsfun, tspan, y0, options, varargin)
    warning on
    try,        h=options.(matchfieldname(options, 'InitialStep'));
    catch,      error(['Option InitialStep must be supplied']);    end
    options.RelTol =1e-3; 
    try,        options.RelTol=options.(matchfieldname(options, 'RelTol'));
    catch,      options.RelTol =1e-3;    end
    % check whether we are going forward or backward in time
    if (tspan(2)<tspan(1)),h=-abs(h);end
    nsteps = abs(ceil ((tspan(2) - tspan(1)) / h));
    ts = zeros(nsteps+1,1);
    ys = zeros(nsteps+1,length(y0));
    t = tspan(1);
    ts(1) = t;
    if (size(y0,2) ~=1), y0=y0';end
    ys(1,:) = y0';
    fac = [];
    options.threshold =1e-6+y0*1e-6;
    vectorized = 0;
    S = [];% generate full Jacobian matrix
    g=[];
    [dFdy,fac,g,nfevals,nfcalls] = numjac(rhsfun,t,y0,feval(rhsfun, t, y0, varargin{:}),options.threshold,fac,vectorized,S,g,varargin{:});
    options.drhsfundy =dFdy;
    for step = 2:nsteps+1
        if (abs(tspan(2)-ts(step-1))<abs(h)),h=tspan(2)-ts(step-1);end
        ts(step) = ts(step-1) + h;
        ys(step,:) = (Trap_Newton (rhsfun, h, ts(step-1), ys(step-1,:)', options, varargin{:}))';
    end
    return;
end

%
% Do one trapezoidal step for the equation dy/dt = f(t,y).
% 
% y_n+1 = y_n + h/2 * (f(t_n+1, y_n+1) +f(t_n, y_n))
%
% Arguments: 
%   rhsfun = rhs function f(t,y)
%   h      = step length
%   tn     = time t_n
%   yn     = function value y_n
%   options= see documentation for odetrap
%   
function [yn1] = Trap_Newton (rhsfun, h, tn, yn, options, varargin)
    tp  = tn + h;
    yp  = yn + h * feval (rhsfun, tp, yn, varargin{:});
    fn=feval(rhsfun, tn, yn, varargin{:});
    yn1 =yp;
    maximum_iteration = 12;
    A = -(h/2) *options.drhsfundy +eye (length(yp));
    for iteration = 1: maximum_iteration
        nyn1 =yn1 - A\(yn1 - yn - (h/2) * (feval(rhsfun, tp, yn1, varargin{:}) + fn));
        yn1increment=nyn1-yn1;
        if norm(yn1increment,inf)<options.RelTol*norm(max(max(abs(yn1)), options.threshold),inf)
            yn1 =nyn1; %iteration, 
            return;
        end
        yn1 =nyn1;
    end
    warning(['odetrap: Not converged with ' num2str(maximum_iteration) ' iterations; Max increment ' num2str(norm(max(max(abs(yn1increment)), options.threshold),inf))]);
end

% This implementation of the trapezoidal step solves for yn1 using
% objective-function minimization, and while being very general, it is also
% quite expensive and slow.
function [yn1] = Trap (rhsfun, h, tn, yn, options, varargin)
    tp  = tn + h;
    yp  = yn + h * feval (rhsfun, tp, yn, varargin{:});
    fn=feval(rhsfun, tn, yn, varargin{:});
    % minimize this function to find yn1
    function v = f(yn1, rhsfun, h, tp, yn, varargin)
        v = yn1 - yn - (h/2) * (feval(rhsfun, tp, yn1, varargin{:}) + fn);
        v = v'*v;
    end
    options =optimset('TolFun',((yn+yp)'*(yn+yp))/1e6,'TolX',h/1e6); % ,'Display','iter'
    yn1 = fminsearch(@f, yp, options, rhsfun, h, tp, yn, varargin{:});
end

