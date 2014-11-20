% Class of a one-point hexahedron quadrature rule.
%
% Parameters:
%      none.
% 
%
classdef hex_4pt_rule
    
    properties
        rule = [];
    end
    
    methods
        
        function self = hex_4pt_rule (Parameters)
            self.rule = gauss_rule(struct('dim',3, 'order',2));
        end
        
        % Get the number of quadrature points of the rule
        function val= npts (self)
            val = length(self.rule.weights)/2;
        end
        
        % Get the weights of quadrature points of the rule
        function val=  weights (self)
            val = self.rule.weights(1:(length(self.rule.weights)/2))*2;
        end
        
        % Get the parametric coordinates of quadrature points of the rule
        function val= param_coords (self)
            val = self.rule.param_coords([1,8,3,6],:);
        end
        
    end
    
end
