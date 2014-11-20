classdef tri_rule
% Class of triangular quadrature rule.
%
% Used for integration on the standard triangle.
    properties
        param_coords = [];% parametric coordinates of the quadrature points
        weights = [];  %  weights of the quadrature points
    end
    
    methods
    
        function self = tri_rule (Parameters)
        % Constructor.
            % Parameters=structure with the mandatory fields
            %      npts=number of points (1-- one-point rule, 3 -- three-point rule,
            %           6 -- six point rule, 10 -- Strang 10 point, order 13, degree of precision 7, rule). 
            % 
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin < 1
                return
            end
            npts = Parameters.npts;
            switch npts
                case 1 % integrates exactly linear polynomials
                    self.param_coords(1,:) = [1/3 1/3];
                    self.weights = [1]/2;
                case 3 % integrates exactly quadratic polynomials
                    self.param_coords = [ 2/3 1/6; 1/6 2/3; 1/6 1/6 ];
                    self.weights = [1/3 1/3 1/3]/2;
                case 6 % integrates exactly quartic polynomials
                    self.param_coords = [ 0.816847572980459 0.091576213509771;...
                        0.091576213509771 0.816847572980459; ...
                        0.091576213509771 0.091576213509771;...
                        0.108103018168070 0.445948490915965;...
                        0.445948490915965 0.108103018168070;...
                        0.445948490915965 0.445948490915965];
                    self.weights = [0.109951743655322*[1, 1, 1] 0.223381589678011*[1, 1, 1] ]/2;
                case 10 % integrates exactly quartic polynomials
                    self.param_coords = [0.333333333333333  0.333333333333333
                        0.479308067841923  0.260345966079038
                        0.260345966079038  0.479308067841923
                        0.260345966079038  0.260345966079038
                        0.869739794195568  0.065130102902216
                        0.065130102902216  0.869739794195568
                        0.065130102902216  0.065130102902216
                        0.638444188569809  0.312865496004875
                        0.638444188569809  0.048690315425316
                        0.312865496004875  0.638444188569809
                        0.312865496004875  0.048690315425316
                        0.048690315425316  0.638444188569809
                        0.048690315425316  0.312865496004875
                        ];
                    self.weights = [ -0.149570044467670
                        0.175615257433204
                        0.175615257433204
                        0.175615257433204
                        0.053347235608839
                        0.053347235608839
                        0.053347235608839
                        0.077113760890257
                        0.077113760890257
                        0.077113760890257
                        0.077113760890257
                        0.077113760890257
                        0.077113760890257
                        ]'/2;
                    
                otherwise
                    error('Unsupported number of points');
            end
        end   
    
        function val=npts (self)
        % Get the number of quadrature points
            val = length(self.weights);
        end
        
    end
    
end
