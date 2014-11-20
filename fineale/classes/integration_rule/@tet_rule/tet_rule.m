% Class of a quadrature rule for tetrahedra.
%
% function self = tet_rule (varargin)
%
% Parameters=structure with the mandatory attributes
%      npts=number of points (1-- one-point rule, 4 -- four-point rule).
%
%
classdef tet_rule
    
    properties
        param_coords = [];
        weights = [];
    end
    
    methods
        
        function self = tet_rule (Parameters)
            if nargin < 1
                return
            end
            npts = Parameters.npts;
            switch npts
                case 1 % integrates exactly linear polynomials
                    self.param_coords = ...
                        [0.25,0.25,0.25];
                    self.weights = [1/6];
                case 4 % integrates exactly quadratic polynomials
                    self.param_coords = ...
                        [[0.13819660,0.13819660,0.13819660];...
                        [0.58541020,0.13819660,0.13819660];...
                        [0.13819660,0.58541020,0.13819660];...
                        [0.13819660,0.13819660,0.58541020]];
                    self.weights = [ 0.041666666666666666667,  0.041666666666666666667,  0.041666666666666666667,  0.041666666666666666667];
                case 5 %   Zienkiewicz #3.
                    a =   1.0 / 6.0;
                    b =   0.25;
                    c =   0.5;
                    d = - 0.8;
                    e =   0.45;
                    self.param_coords = ...
                        [[b,b,b];...
                        [c,a,a];...
                        [a,c,a];...
                        [a,a,c];...
                        [a,a,a]];
                    self.weights = [d, e, e, e, e]/6;
                otherwise
                    error('Unsupported number of points');
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
