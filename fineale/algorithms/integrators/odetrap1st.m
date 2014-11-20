function [ts, ys] = odetrap1st (rhsfun, tspan, y0, options, varargin)
% Trapezoidal integrator (first-order form).
  %
% Function to integrate a system of first-order 
% initial-value ODE's by the trapezoidal method.
% 
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
    h = options.InitialStep(1);
    K=options.K;
    M=options.M;
    F0=zeros(size(M,1));
    F=@(t)F0;
    if isfield( options,'F')
        F=options.F;
    end
    maxit=12;
    if isfield( options,'maxit')
        maxit=options.maxit;
    end
    nsteps = ceil ((tspan(2) - tspan(1)) / h);
    ts = zeros(nsteps+1,1);
    ys = zeros(nsteps+1,length(y0));
    n=length(y0)/2;
    t = tspan(1);
    ts(1) = t;
    ys(1,:) = y0';
    F0 =F(ts(1));
    A = [zeros(n),eye(n);-M\K,zeros(n)];
    I = eye(length(y0));
    for step = 2:nsteps+1
        if (ts(step-1) + h >tspan(2)), h=tspan(2)-ts(step-1); end
        F1 =F(ts(step-1) + h);
        y0 =ys(step-1,:)';
        L0=[zeros(n,1);M\F0];
        L1=[zeros(n,1);M\F1];
        y1 =(I-h/2*A)\((I+h/2*A)*y0+h/2*L0+h/2*L1);
        ts(step) = ts(step-1) + h;
        ys(step,:) = y1';
        F0 =F1;
    end
    return;
end
