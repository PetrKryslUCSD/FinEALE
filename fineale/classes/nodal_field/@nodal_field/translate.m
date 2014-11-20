function self = translate(self, c)
% Apply translation to each nodal parameter.
%
% function self = translate(self, c)
%
    [nfens,dim] = size(self.values);
    for i=1:nfens
        self.values(i,:)=self.values(i,:)+c;
    end
end
