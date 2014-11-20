function val = gather_by_dofnums (self, dofnums)
% Gather values from the field corresponding to the degree-of-freedom numbers.
%
% function val = gather_by_dofnums (self, dofnums)
%
% Gather values from the field corresponding
% to the degree-of-freedom numbers.
%
    n = length(dofnums);
    val = zeros(n,1);
    [nfens,dim] = size(self.values);
    for i=1:nfens
        for j=1:dim
            en = self.dofnums(i,j);
            for k=1:n
                if (dofnums(k) == en)
                    val(k) = self.values(i,j);
                end
            end
        end
    end
    return;
end

