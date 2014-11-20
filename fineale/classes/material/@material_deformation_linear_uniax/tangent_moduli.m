% Calculate the material stiffness matrix.
%
% function D = tangent_moduli(self, context)
%
%   Call as:
%     D = tangent_moduli(m, context)
%   where
%     m=material
%    context=structure with optional field:
%
% the output arguments are
%     D=matrix 1x1
%
function D = tangent_moduli(self, context)
    D = get(self.property,'E');
    return;
end

