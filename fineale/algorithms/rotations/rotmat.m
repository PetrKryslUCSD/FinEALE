% Compute rotation matrix from a rotation vector or from the associated
% skew matrix.
%
% function R = rotmat(theta)
%
% Compute rotation matrix from a rotation vector or from the associated
% skew matrix theta.
function R = rotmat(theta)
[m,n]=size(theta);
if m == 3 && n == 3
    if (norm(theta' + theta) > 1e-9)
        error('theta must be a skew-symmetric matrix');
    end
    thetatilde=theta;
else
    thetatilde=[  0     -theta(3) theta(2); ...
                theta(3)    0    -theta(1); ...
               -theta(2) theta(1)   0];
end
%R = expm(thetatilde);
a=[ -thetatilde(2,3); thetatilde(1,3); -thetatilde(1,2)];
na = norm(a);
if (na == 0)
    R = eye(3,3);
else
    a=a/na;
    ca = cos(na);
    sa = sin(na);
    R = ca * (eye(3,3)-(a*a')) + sa/na*thetatilde + (a*a');
end
return;

