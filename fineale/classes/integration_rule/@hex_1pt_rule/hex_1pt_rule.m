% Class of a one-point hexahedron quadrature rule.
%
% Parameters:
%      none.
% 
%
classdef hex_1pt_rule
    
    properties
        rule = [];
    end
    
    methods
        
        function self = hex_1pt_rule (Parameters)
            self.rule = gauss_rule(struct('dim',3, 'order',1));
        end
        
        % Get the number of quadrature points of the rule
        function val= npts (self)
            val = length(self.rule.weights);
        end
        
        % Get the weights of quadrature points of the rule
        function val= weights (self)
            val = self.rule.weights;
        end
        
        % Get the parametric coordinates of quadrature points of the rule
        function val= param_coords (self)
            val = self.rule.param_coords;
        end
        
    end
    
end
