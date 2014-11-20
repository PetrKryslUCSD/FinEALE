function vec = gather_sysvec (self)
% Gather values from the field for the whole system vector.
%
% function vec = gather_sysvec (self)
%
%   Call as:
%       v=gather_sysvec(f)
%
    [nfens,dim] = size(self.values);
    vec=zeros(self.nfreedofs,1);
    for i=1:nfens
        for j=1:dim
            en = self.dofnums(i,j);
            if (en > 0 & en <= self.nfreedofs)
                vec(en)=self.values(i,j);
            end
        end
    end
    return;
end

