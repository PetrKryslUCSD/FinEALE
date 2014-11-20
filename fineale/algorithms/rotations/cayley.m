% Compute the Cayley transform.
%
% function C = cayley(theta)
%
% theta = rotation parameters
function C = cayley(theta)
[m,n]=size(theta);
if m == 3  &&  n == 3
    if (norm(theta' + theta) > 1e-9)
        error('theta must be a skew-symmetric matrix');
    end
    thetatilde=theta;
else
    thetatilde=[  0     -theta(3) theta(2); ...
                theta(3)    0    -theta(1); ...
               -theta(2) theta(1)   0];
end
C=(eye(3,3)+1/2*thetatilde)*inv(eye(3,3)-1/2*thetatilde);
return;

