classdef femm_deformation_linear < femm_deformation
    % Class for the small displacement, small strain deformation model.
    %
    
    properties
    end
    
    
    methods % constructor
        
        function self = femm_deformation_linear (Parameters)
            % Constructor.
            % Parameters: those recognized by femm_deformation 
            %
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin <1
                Parameters = struct([]);
            end
            self = self@femm_deformation(Parameters);
            if (nargin <1)  || isempty(Parameters)
                return
            end
        end
        
    end
    
end

