function retobj = clone(self,name)
% Make a new nodal field by cloning another field
%
% function retobj = clone(fld,name)
%
%   Call as:
%     clonedf = clone(f, name)
%   where
%     f=field to clone
%     name=name of the new field
%
    retobj =  self;
    retobj.name = name;
    return;
end

