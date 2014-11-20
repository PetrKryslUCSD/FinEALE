% Update material state. 
%
% function [out, newms] = update (self, ms, context)
%
% Update material state.  Return the updated material state, and the
% requested quantity (default is the stress).
%   Call as:
%     [out,newms] = update(m, ms, context)
%  where
%     m=material
%     ms = material state
%     context=structure
%        with mandatory fields
%           gradU= displacement gradient
%        and optional fields
%           output=type of quantity to output, and interpreted by the
%           particular material; [] is returned when the material does not
%           recognize the requested quantity to indicate uninitialized
%           value.  It can be tested with isempty ().
%              output ='Cauchy' - Cauchy stress; this is the default
%                      when output type is not specified.
%              output ='2ndPK' - 2nd Piola-Kirchhoff stress;
%           It is assumed that stress is output in 6-component vector form.
%   The output arguments are
%     out=requested quantity
%     newms=new material state; don't forget that if the update is final
%           the material state newms must be assigned and stored.  Otherwise
%           the material update is lost!
%
function [out, newms] = update (self, ms, context)
    strain = context.strain;
    D  = tangent_moduli (self, context);
    tSigma = thermal_stress(self,context);
    stress = D * Ev + tSigma;
    if isfield(context,'output')
        switch context.output
            case 'Cauchy'
                out = stress;
            otherwise
                out = [];
        end
    else
        out = stress;
    end
    newms = ms;
    return;
end
