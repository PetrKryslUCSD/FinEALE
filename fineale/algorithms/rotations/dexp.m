% Compute the derivative of the exponential map operator.
%
% function m=dexp(x)
%
% x= rotation vector
function m=dexp(x)
phi=sqrt(x'*x);
if (phi == 0)
    m=eye(3,3);
else
    xs=skewmat(x);
    c1=(1-cos(phi))/phi/phi;
    c2=(1-sin(phi)/phi)/phi/phi;
    m=eye(3,3) + c1*xs + c2*(xs*xs);
end
