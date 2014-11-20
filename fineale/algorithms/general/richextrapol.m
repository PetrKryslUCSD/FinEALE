function [xestim, beta, c, residual] = richextrapol(xs,hs)
% Richardson extrapolation.
%
% function [xestim, beta, c, residual] = richextrapol(xs,hs)
%
% Richardson extrapolation. This function is applicable only to fixed ratio
% between the mesh sizes, hs(1)/hs(2) = hs(2)/hs(3).
% xs = the calculated quantities,
% hs = the mesh sizes
%
% Returns
% xestim= estimate of the asymptotic solution from the data points in the xs array
% beta= convergence rate
% c = constant in the estimate "error=c*h^beta"
% residual = residual after equations from which the above quantities were
%      solved (this is a measure of how accurately was the system solved)
    if abs(hs(1)/hs(2) - hs(2)/hs(3)) > 1e-6
        warning ([' Ratio between mesh sizes should be fixed: hs(1)/hs(2) - hs(2)/hs(3) = 0; I got ' num2str(hs(1)/hs(2) - hs(2)/hs(3))]);
    end
    nxs=xs./xs(1);
    try
    xestim=eval(vpa(-vpa(nxs(1)*nxs(3)-nxs(2)^2)/vpa(2*nxs(2)-nxs(1)-nxs(3))))*xs(1);
    catch
    xestim=((-(nxs(1)*nxs(3)-nxs(2)^2)/(2*nxs(2)-nxs(1)-nxs(3))))*xs(1);
    end
    if (xestim-xs(1)) <= 0
        beta=log((xestim-xs(2))./(xestim-xs(3)))/log(hs(2)/hs(3));
    else
        beta=log((xestim-xs(1))./(xestim-xs(3)))/log(hs(1)/hs(3));
    end
    % just to check things, calculate the residual
    c=(xestim-xs(1))/hs(1)^beta;
    for I =1:3
        residual(I) =(xestim-xs(I))-c*hs(I)^beta;% this should be close to zero
    end
    return;
end
