classdef trapezoidal_rule
% Class of trapezoidal quadrature rule.
%
% The rule is applicable for a tensor product of  intervals -1 <=x<= +1.
% 


    properties
        dim = [];%  number of space dimensions (1,2,3)
        param_coords = [];% parametric coordinates of the quadrature points
        weights = [];  %  weights of the quadrature points
    end
    
    methods
    
        function self = trapezoidal_rule (Parameters)
        % Constructor.
            % Parameters: 
            %   dim=number of  space dimensions
            %            3=3-D rule, 2=2-D  rule, 1=1-D rule (the default).
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin < 1
                return
            end
            self.dim = Parameters.dim;
            switch self.dim
                case 1
                    self.param_coords = [-1; 1];
                    self.weights = [1, 1]';
                case 2
                    self.param_coords = [-1,-1; 1,-1; 1, 1; -1, 1];
                    self.weights = [1, 1, 1, 1]';
                case 3
                    self.param_coords = [-1,-1,-1; 1,-1,-1; 1, 1,-1;-1, 1,-1;-1,-1, 1; 1,-1, 1; 1, 1, 1;-1, 1, 1;];
                    self.weights = [1, 1, 1, 1, 1, 1, 1, 1]';
                otherwise
                    error('Unsupported Dimension');
            end
        end
    
        function val= npts (self)
        % Get the number of quadrature points of the rule
            val = length(self.weights);
        end
        
        function val= get.weights (self)
        % Get the weights of quadrature points of the rule
            val = self.weights;
        end
        
        function val= get.param_coords (self)
        % Get the parametric coordinates of quadrature points of the rule
            val = self.param_coords;
        end
        
    end
    
end
 