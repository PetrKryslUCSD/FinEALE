% Compute a skew matrix from its axial vector Theta.
%
% function S = skewmat(theta)
%
function S = skewmat(theta)
    [m,n]=size(theta);
    if m*n == 3
        S=[  0     -theta(3) theta(2); ...
            theta(3)    0    -theta(1); ...
            -theta(2) theta(1)   0];
    elseif  m*n == 2
        S=[-theta(2) theta(1)];
    else
        error('theta must be a 3-vector');
    end
    return;
end
