function self = scatter_sysvec (self, vec)
% Scatter values to the field from a system vector.
%
% function retobj = scatter_sysvec (self, vec)
%
%   Call as:
%     f=scatter_sysvec(f, vec)
%   where
%     f=field (note that we have to assign to f: f changes inside this function)
%     vec=system vector (an array of doubles), indexed by free degree of freedom numbers
%
    [nfens,dim] = size(self.values);
    for i=1:nfens
        for j=1:dim
            dn = self.dofnums(i,j);
            if (dn > 0) && (dn <= self.nfreedofs)
                self.values(i,j) = vec(dn);
            end
        end
    end
end
