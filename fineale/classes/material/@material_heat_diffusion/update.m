function [out, newms] = update (self, ms, context)
% _Update_ material state.  
%
% function [out, newms] = update (self, ms, context)
%
% Update material state.  Return the updated material state, and the
% requested quantity (default is the heat flux).
%    
% Arguments
%     m=material
%     ms = material state
%     context=structure
%       with mandatory fields
%         gradtheta= temperature gradient in the local material
%           directions (which may be the same as the global coordinate
%           directions)
%         Rm = material orientation  matrix, may be supplied as empty  if it 
%           corresponds to an identity.   
%       and optional fields
%         output=type of quantity to output, and interpreted by the
%           particular material; [] is returned when the material does not
%           recognize the requested quantity to indicate uninitialized
%           value.  It can be tested with isempty ().
%              output ='flux' - heat flux vector in the local material 
%                      directions; this is the default
%                      when output type is not specified.
%         outputRm = orientation matrix  of the coordinate system in which 
%           the output should be calculated;  if this matrix is not supplied, 
%           it is assumed that the output is to be provided  in the local 
%           material coordinate system; otherwise the output is first 
%           transformed to the global coordinate system, and then to the 
%           output coordinate system.
% Output
%     out=requested quantity
%     newms=new material state; don't forget that if the update is final
%           the material state newms must be assigned and stored.  Otherwise
%           the material update is lost!
%
    kappa = self.property.thermal_conductivity;
    gradtheta = reshape (context.gradtheta, [],1);% make column vector
    flux = - kappa * gradtheta;% in local material orientation coordinates
    if (isfield(context,'outputRm'))%  output coordinate system supplied?
        if (isfield(context,'Rm')) && ( ~isempty(context.Rm) )
            flux = (context.Rm*flux);% in global coordinate system
        end
        flux = context.outputRm'*flux;% in output coordinate system
    end    
    if isfield(context,'output')
        switch context.output
            case 'flux'
                out = flux;
            otherwise
                out = [];
        end
    else
        out = flux;
    end
    newms = ms;
    return;
end
