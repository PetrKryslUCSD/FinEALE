% Class of a special six-point hexahedron quadrature rule.
%
% Parameters:
%      none.
% 
% The integration scheme is based on the paper
% Sauer G (1999) Alternative reduced integration avoiding spurious modes for 8-node quadrilateral and 20-node hexahedron finite elements. Forschung Im Ingenieurwesen-Engineering Research 65: 131-135.
%
classdef hex_24pt_rule
    
    properties
        param_coords = [];
        weights = [];
    end
    
    methods
        
        function self = hex_24pt_rule (Parameters)
            if (nargin<1)
                return
            end
            a= Parameters.a;
            b= Parameters.b;
            %
            C = [...
            [-a,-a,-b];...
            [-a,-a,+b];...
            [+a,-a,-b];...
            [+a,-a,+b];...
            [+a,+a,-b];...
            [+a,+a,+b];...
            [-a,+a,-b];...
            [-a,+a,+b];...
            ];
            self.param_coords= [C(:,[1,2,3]); C(:,[1,3,2]); C(:,[3,1,2])];;
            self.weights = [ 8/24*ones(1,24)]';
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
