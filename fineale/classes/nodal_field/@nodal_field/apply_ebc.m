function self = apply_ebc (self)
% Apply EBCs (essential boundary conditions).
%
% function retobj = apply_ebc (self)
%
% This function applies existing boundary conditions.  In other
% words, it sets value=fixed_values for all parameters and components
% of parameters that are being fixed as essential boundary
% conditions (EBCs).
%
% To set the prescribed values and to mark components as being prescribed,
% use set_ebc().
%
% Don't forget to assign the field to itself, because this function
% changes the internal state of self.
%    Call as:
%      f = apply_ebc(f);
%
    [nfens,dim] = size(self.values);
    for i=1:nfens
        for j=1:dim
            if (self.is_fixed(i,j))
                self.values(i,j) = self.fixed_values(i,j);
            end
        end
    end
end
