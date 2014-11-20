function self = scatter_by_dofnums (self, dofnums, val)
% Scatter values to the field by  degree-of-freedom numbers.
%
% function retobj = scatter_by_dofnums (self, dofnums, val)
%
% Scatter values to the field corresponding
% to the degree of freedom numbers.
%   Call as:
%     f=scatter_by_dofnums(f, dofnums, val)
%   where
%     f=field (note that we have to assign to f: f changes inside this function)
%     dofnums=array of degree-of-freedom numbers
%     val=array of values to set (i.e. double numbers)
%
    n = length(dofnums);
    [nfens,dim] = size(self.values);
    for i=1:nfens
        for j=1:dim
            en = self.dofnums(i,j);
            for k=1:n
                if ((dofnums(k) > 0) & (dofnums(k) == en))
                    self.values(i,j) = val(k);
                end
            end
        end
    end
end
