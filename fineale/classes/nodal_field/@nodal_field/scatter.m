function self = scatter (self, fenids, val)
% Distribute (scatter) values to the field by the finite element node indices.
%
% function retobj = scatter (self, fenids, val)
%
% Scatter values to the field corresponding
% to the finite element node indices.
%   Call as:
%     f=scatter(f, fenids, val)
%   where
%     f=field (note that we have to assign to f: f changes inside this function)
%     fenids=array of node identifiers
%     val=an array of doubles, size(val) == [length(fenids),dim]
%
    [nfens,dim] = size(self.values);
    if (isempty(fenids))
        fenids = (1:nfens);
    end
    n = length(fenids);
    for i=1:n
        k = fenids(i);
        self.values(k,:) = val(i,:);
    end
end

