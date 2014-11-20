% Class of a special six-point hexahedron quadrature rule.
%
% Parameters:
%      none.
%
classdef hex_6pt_rule
    
    properties
        dim = [];
        order = [];
        param_coords = [];
        weights = [];
    end
    
    methods
        
        function self = hex_6pt_rule (Parameters)
            self.param_coords = ...
            [[-1, 0, 0];...
            [1, 0, 0];...
            [0,-1,0];...
            [0,1,0];...
            [0, 0,-1];...
            [0, 0, 1]];
            self.weights =4/3*ones(6, 1);
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
