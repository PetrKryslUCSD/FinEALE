function self = slice(self,idc,name)
% Make a new nodal_field by taking a slice of another field.
%
% function self = slice(fld,idc,name)
%
%   Call as:
%     nf = slice(f, idc, name)
%   where
%     f=field to clone
%     idc = array of indexes of components; the new field consists of the
%     specified components taken out of the nodal parameters of the
%     original field. All attributes are extracted from the original field:
%     subset of the values, equation numbers,is_fixed attribute, and
%     fixed_values attribute.
%     name=name of the new field
%
    self.name=name;
    dim  = length(idc);
    self.values =self.values(:,idc);
    self.dofnums =self.dofnums(:,idc);
    self.is_fixed = self.is_fixed(:,idc);
    self.fixed_values = self.fixed_values(:,idc);
    self.nfreedofs = self.nfreedofs;
end

