% Class of the 'point' quadrature rule (used exclusively for zero-dimensional manifolds).
%
% Parameters; none
% Example:
%      ir = point_rule

classdef point_rule

    properties
        param_coords = 0;
        weights = 1;  
    end
    
    methods
    
        function self = point_rule (Parameters)
            if nargin < 1
                return
            end
        end
    
        % Get the number of quadrature points of the rule
        function val= npts (self)
            val = length(self.weights);
        end
        
        % Get the weights of quadrature points of the rule
        function val= get.weights (self)
            val = self.weights;
        end
        
        % Get the parametric coordinates of quadrature points of the rule
        function val= get.param_coords (self)
            val = self.param_coords;
        end
        
    end
    
end
