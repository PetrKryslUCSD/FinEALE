function self = rotate(self, R)
% Apply rotation matrix to each nodal parameter.
%
% function retobj = rotate(self, R)
%
% Note: makes sense only for 3D geometry fields!
%
    [nfens,dim] = size(self.values);
    if (dim ~= 3)
        error('Dimensiom mismatch! Need 3.');
    end
    for i=1:nfens
        self.values(i,:)=(R*self.values(i,:)');
    end
    return;
end
