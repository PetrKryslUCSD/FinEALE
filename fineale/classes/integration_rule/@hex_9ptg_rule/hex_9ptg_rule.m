% Class of a special 9-point hexahedron quadrature rule.
%
% Parameters:
%      none.
% 
% The integration scheme is based on the paper
% Sauer G (1999) Alternative reduced integration avoiding spurious modes 
% for 8-node quadrilateral and 20-node hexahedron finite elements. 
% Forschung Im Ingenieurwesen-Engineering Research 65: 131-135.
%
classdef hex_9ptg_rule
    
    properties
         param_coords = [];
        weights = [];
    end
    
    methods
        
        function self = hex_9ptg_rule (Parameters)
            rh= sqrt(1/3);
        %
        self.param_coords = [...
        [-rh,-rh,-rh];...
        [+rh,-rh,-rh];...
        [+rh,+rh,-rh];...
        [-rh,+rh,-rh];...
        [-rh,-rh,+rh];...
        [+rh,-rh,+rh];...
        [+rh,+rh,+rh];...
        [-rh,+rh,+rh];...
        [0,0,0];...
        ];
        self.weights = [695/700*ones(1,8),(8-8*695/700)]';
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
