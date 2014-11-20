function fxi = simplex_grid_interpolation(x,fx,tri,xi)
% Interpolation on a simplex grid.
%
% function fxi = simplex_grid_interpolation(x,fx,tri,xi)
%
%   x = dimension m-by-n, representing m points in n-D space (n=2,3). 
%   fx = dimension m-by-1, representing m values of the function F(X).  
%   xi = a vector of size p-by-n, representing p points in the n-D space

    % Find the nearest triangle (t)
    [t,p] = tsearchn(x,tri,xi);

    fxi = NaN*zeros(size(xi,1),1);

    for i = 1:length(fxi)
        if ~isnan(t(i))
            fxi(i) = p(i,:)*fx(tri(t(i),:));
        end
    end
end
