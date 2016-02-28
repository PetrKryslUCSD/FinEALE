% Class of simplex nodal-quadrature rule.
% 
% function retobj = simplex_nq_rule (Parameters)
% 
%   Parameters=structure with fields
%      dim=number of  space dimensions
%            3=3-D rule (four-node tetrahedron), 2=2-D  rule (three-node triangle).
% 
%
classdef simplex_nq_rule
    
    properties
        param_coords = [];
        weights = [];
    end
    
    methods
        
        function self = simplex_nq_rule (Parameters)
            if nargin < 1
                return
            end
            dim = Parameters.dim;
            if (isfield( Parameters, 'order' ))
                order =  Parameters.order;
            else
                order = 1;
            end
            switch dim
                case 1
                    self.param_coords = [-1;+1];
                    self.weights = [1,1];
                case 2
                    self.param_coords = [ 0, 0; 1, 0; 0, 1];
                    self.weights = [1/3 1/3 1/3]/2;
                case 3
                    if (order == 1)
                        self.param_coords = [0, 0, 0; 1, 0, 0; 0, 1, 0; 0, 0, 1];
                        self.weights = [1, 1, 1, 1]/6/ 4;
                    elseif (order == 2)
                        b = 0; %  Rule with zero weight for the mid-sides
                        % b = (1)/6; %  Rule with zero weight for the corners
                        % b = (1.2)/6; % NCC2, the Newton Cotes Closed Rule,
                        % order 10, degree of precision 2.
                        % http://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tet/quadrature_rules_tet.html
                        % Peter Silvester,
                        % Symmetric Quadrature Formulae for Simplexes,
                        % Mathematics of Computation,
                        % Volume 24, Number 109, January 1970, pages 95-100.
                        a = (1- 6*b)/4;
                        self.param_coords = [0     0     0
                                             1     0     0
                                             0     1     0
                                             0     0     1
                                            0.5000         0         0
                                            0.5000    0.5000         0
                                            0    0.5000         0
                                            0         0    0.5000
                                            0.5000         0    0.5000
                                            0    0.5000    0.5000];
                        self.weights = [a*ones(1,4), b*ones(1,6)]/6;
                    else    
                        error('Unsupported order');
                    end
                otherwise
                    error('Unsupported dimension');
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
