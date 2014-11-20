% Calculate vector of thermal stress components.
%
% function v = thermal_stress(self,context)
%
%   Call as:
%     v = thermal_stress(m,context)
%  where
%     m=material
%     context=structure; see the update() method
%
function v = thermal_stress(self,context)
    alphas = get(self.property,'alphas');
    D = tangent_moduli(self, context);
    v = -D*context.dT*alphas(1);
end
