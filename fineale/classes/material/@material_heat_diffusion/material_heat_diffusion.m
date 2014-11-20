classdef material_heat_diffusion < material_base
% Class for heat diffusion models for materials.
%
% This class represents linear diffusion models in materials.
%

    properties
        property= [];% heat diffusion property object
    end
    
    methods
    
        function self = material_heat_diffusion (Parameters)
        % Constructor.
        % Parameters:
        %     property=heat-diffusion property object
        % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin <1
                return
            end
            self.property = Parameters.property;
        end
        
    end
    
end
