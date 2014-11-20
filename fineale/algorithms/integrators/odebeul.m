% Backward Euler integrator.
%
% Function to integrate a system of first-order 
% initial-value ODE's by the backward Euler method.
% 
% One Backward Euler step for the equation dy/dt = f(t,y) reads
% 
% y_n+1 = y_n + h * f(t_n+1, y_n+1)
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
function [ts, ys] = odebeul (rhsfun, tspan, y0, options, varargin)
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
        ys(step,:) = (BEuler_Newton(rhsfun, h, ts(step-1), ys(step-1,:)', options, varargin{:}))';
    end
    return;
end

%
% Do one Backward Euler step for the equation dy/dt = f(t,y).
% 
% y_n+1 = y_n + h * f(t_n+1, y_n+1)
%
% Arguments: 
%   rhsfun = rhs function f(t,y)
%   h      = step length
%   tn     = time t_n
%   yn     = function value y_n
%   
function [yn1] = BEuler_Newton (rhsfun, h, tn, yn, options, varargin)
    tp  = tn + h;
    yp  = yn + h * feval (rhsfun, tp, yn, varargin{:});
    yn1 =yp;
    maximum_iteration = 20;
    A = -(h) *options.drhsfundy +eye (length(yp));
    for iteration = 1: maximum_iteration
        nyn1 =yn1 - A\(yn1 - yn - (h) * feval(rhsfun, tp, yn1, varargin{:}));
        yn1increment=nyn1-yn1; 
        if norm(yn1increment,inf)<options.RelTol*norm(max(max(abs(yn1)), options.threshold),inf)
            yn1 =nyn1; return;
        end
        yn1 =nyn1;
    end
    yn1increment(isnan(yn1increment)) = realmax;
    warning(['odebeul: Not converged with ' num2str(maximum_iteration) ' iterations; Max increment ' num2str(norm(max(max(abs(yn1increment)), options.threshold),inf))]);
end

% function [yn1] = BEuler (rhsfun, h, tn, yn, varargin)
%     yn1 = BEuler_fminsearch (rhsfun, h, tn, yn, varargin{:});
% end
% 
% function [yn1] = BEuler_fminsearch (rhsfun, h, tn, yn, varargin)
%     tp  = tn + h;
%     yp  = yn + h * feval (rhsfun, tp, yn, varargin{:});
%     options =optimset('TolFun',((yn+yp)'*(yn+yp))/1e16,'MaxIter',100); % ,'Display','iter','TolX',0
%     [yn1,fval] = fminsearch(@beulerf, yp, options, rhsfun, h, tp, yn, varargin{:}); % optimset('disp','iter')
% end
% 
% function v = beulerf(yn1, rhsfun, h, tp, yn, varargin)
%     v = yn1 - yn - h * feval (rhsfun, tp, yn1, varargin{:});
%     v = v'*v ;
% end
% 
% function [yn1] = BEuler_fixed_point (rhsfun, h, tn, yn, varargin)
%     tol = 1e-6;
%     yp  = yn; % initial guess
%     tp  = tn + h;
%     it  = 1;
%     while (1)
%         prevyp = yp;
%         yp = yn + h * feval (rhsfun, tp, yp, varargin{:});
%         if (norm(yp - prevyp) < tol)
%             % fprintf (1, '# it = %d\n', it);
%             break;
%         end;
%         it = it + 1;
%     end
%     yn1 = yp;
% end