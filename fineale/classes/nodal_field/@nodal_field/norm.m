function n = norm(self,varargin)
% Compute the norm of the field values array.
%
% function n = norm(self,varargin)
%
%       The array is treated as a vector.  In other words,
%       if you gather a system vector equation by equation and
%       compute its norm, it is what you'd get with norm(field).
%
    which_norm = 2;
    if (nargin >= 2)
        which_norm = varargin{1};
    end
    switch which_norm
        case 2
            n = 0;
            for i = 1:self.dim
                n = n + self.values(:,i)' * self.values(:,i);
            end
            n = sqrt(n);
        case 1
             n = norm(self.values,1);
        case inf
             n = norm(self.values,inf);
        otherwise
            error(['Unsupported norm requested: ' num2str(which_norm) '!']);
    end
    return;
end


