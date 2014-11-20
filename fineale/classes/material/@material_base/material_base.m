classdef material_base
% Class of an abstract material.  Base class for all other materials. 
%
%

    properties
    
    end
    
    methods
    
        function self = material_base(Parameters)
        % Constructor.
        % Parameters:  none.
        %
        % This class represents materials. It is abstract.
        % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if (nargin < 1)
                Parameters =struct( [] );
            end
        end
        
        function v = newmatstate (self)
        % Create a new material state.
        %
        % function v = newmatstate (self)
        %
        %     self=material descendent
        %
        % Create data structure to hold the material state at the quadrature point.
        %
        % No implementation in the base class.
            error ('Not implemented in the abstract class');
        end
        
        function [out, newms] = update (self, ms, context) 
        % _Update_ material state. 
        %
        % function [out, newms] = update (self, ms, context)
        %
        %     self=material descendent
        %     ms = material state
        %     context=structure
        %
        % Update material state.  Return the updated material state, and the
        % requested quantity (default is the heat flux).
        %   
        % No implementation in the base class.
            error ('Not implemented in the abstract class');
        end

        function draw (self, gv, context)    
        % Produce a graphic representation of data at the integration point.
        %
        % function draw (self, gv, context)
        %
        %
        % Input arguments
        % self = self
        % gv = graphic viewer
        % context = struct with  appropriate context data
        %  
        % No implementation in the base class.
            error ('Not implemented in the abstract class');
        end
    
    end
    
end
