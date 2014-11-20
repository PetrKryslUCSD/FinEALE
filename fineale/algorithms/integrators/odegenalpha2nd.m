% Generalized-alpha method for a second order system.
%
% Function to integrate a system of second-order (but cast as a first order
% system) initial-value ODE's by the generalized-alpha method
% Title: A TIME INTEGRATION ALGORITHM FOR STRUCTURAL DYNAMICS WITH IMPROVED NUMERICAL DISSIPATION - THE GENERALIZED-ALPHA METHOD
% Author(s): CHUNG J, HULBERT GM
% Source: JOURNAL OF APPLIED MECHANICS-TRANSACTIONS OF THE ASME   Volume: 60   Issue: 2   Pages: 371-375   Published: JUN 1993
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
function [ts, ys] = odegenalpha2nd (rhsfun, tspan, y0, options, varargin)
    h = options.InitialStep(1);
    if ~isfield( options,'RelTol')
        options.RelTol =1e-3;
    else
        if isempty(options.RelTol)
            options.RelTol =1e-3;
        end
    end
    K=options.K;
    M=options.M;
    C=zeros(size(M));
    if isfield( options,'C')
        C=options.C;
    end
    F0=zeros(size(M,1));
    F=@(t)F0;
    if isfield( options,'F')
        F=options.F;
    end
    % HHT-alpha: alphm = 0, alphf = (1-rhoinf)/(1+rhoinf), rhoinf~1.0--0.5
    % WBZ-alpha: alphm = -(1-rhoinf)/(1+rhoinf), alphf = 0, rhoinf~1.0--0.5
    % Generalized: alphf = (alphm+1)/3
    alphm=0;
    if isfield( options,'alphm')
        alphm=options.alphm;
    end
    alphf=0;
    if isfield( options,'alphf')
        alphf=options.alphf;
    end
    gam= 1/2-(alphm-alphf);
    if isfield( options,'gam')
        gam=options.gam;
    end
    %1/4*(1/2+Newmarkg)^2
    bet = 1/4*(1/2+gam)^2 -1/2*(alphm-alphf);
    if isfield( options,'bet')
        bet=options.bet;
    end
    maxit=1;
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
    % dynamic stiffness
    S=((1-alphm)*M+h*gam*(1-alphf)*C+h^2*bet*(1-alphf)*K);
    % initial setup
    tn=ts(1);
    dn=ys(1,1:n)';
    vn=ys(1,n+1:end)';
    an=M\(F(tn)-K*dn-C*vn);
    for step = 2:nsteps+1
        if (ts(step-1) + h >tspan(2)), h=tspan(2)-ts(step-1); end
        tn=ts(step-1);
        dn=ys(step-1,1:n)';
        vn=ys(step-1,n+1:end)';
        tn1=tn + h;
        %         M*an1am+C*vn1af+K*dn1af=Fn1af
        %         dn1af=(1-alphf)*dn1+alphf*dn;
        %        vn1af=(1-alphf)*vn1+alphf*vn;
        %         an1am= (1-alphm)*an1+alphm*an;
        %       M*(1-alphm)*an1+M*alphm*an+C*((1-alphf)*(vn+h*((1-gam)*an+gam*an1))+alphf*vn)+K*((1-alphf)*(dn+h*vn+h^2*((1/2-bet)*an+bet*an1))+alphf*dn)
        Fn1af=F((1-alphf)*tn1+alphf*tn);
        R=-(M*alphm*an+C*((1-alphf)*(vn+h*((1-gam)*an))+alphf*vn)+K*((1-alphf)*(dn+h*vn+h^2*((1/2-bet)*an))+alphf*dn))+Fn1af;
        an1=S\R;
        dn1 =dn+h*vn+h^2*((1/2-bet)*an+bet*an1);
        vn1=vn+h*((1-gam)*an+gam*an1);
        an =an1;
        ts(step) = tn1;
        ys(step,:) = [dn1',vn1'];
    end
    return;
end
