% Compute the inverse of the derivative of the exponential map.
%
% function m=dexpinv(x)
%
function m=dexpinv(x)
phi=sqrt(x'*x);
if (phi == 0)
    m=eye(3,3);
else
    xs=skewmat(x);
    c1=-1/2;
    c2=-((phi/2)*cot(phi/2)-1)/phi/phi;
    m=eye(3,3) + c1*xs + c2*(xs*xs);
end
