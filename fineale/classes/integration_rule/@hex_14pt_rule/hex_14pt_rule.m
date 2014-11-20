% Class of a special six-point hexahedron quadrature rule.
%
% Parameters:
%      none.
%
% The integration scheme is based on the paper
% Hoit M, Krishnamurthy K (1995)
% A 14-POINT REDUCED INTEGRATION SCHEME FOR SOLID ELEMENTS.
% Computers & Structures 54: 725-730.
%
classdef hex_14pt_rule
    
    properties
        % syms Wa real
        % syms Wb real
        % syms a real
        % syms b real
        % ss=solve(8*Wa+6*Wb-8,8*a^2*Wa+2*b^2*Wb-8/3,8*a^4*Wa+2*b^4*Wb-8/5,8*a^4*Wa-8/9)
        % ss.a, ss.Wa, ss.b, ss.Wb
        %
        % Exact solution of the above equations
        a=627^(1/2)/33;
        Wa=121/361;
        b=570^(1/2)/30;
        Wb=320/361;
        %
        % syms Wb real
        % syms a real
        % syms b real
        % Wa = 0.9999
        % ss=solve(8*Wa+6*Wb-8,8*a^2*Wa+2*b^2*Wb-8/3,8*a^4*Wa-8/9)
        % ss.a,  ss.b, ss.Wb
        %
        % Special solution for Wa = 0.9999
        % a=(3333^(1/2)*(-((10000 - 300*1111^(1/2))^(1/2) - 100)*((10000 - 300*1111^(1/2))^(1/2) + 100))^(1/2))/9999;
        % Wa=0.9999;
        % b=(10000 - 300*1111^(1/2))^(1/2);
        % Wb=1/7500;
        %
        % syms Wa real
        % syms Wb real
        % syms a real
        % syms b real
        % Wa = 0.9
        % ss=solve(8*Wa+6*Wb-8,8*a^2*Wa+2*b^2*Wb-8/3,8*a^4*Wa-8/9)
        % ss.a,  ss.b, ss.Wb
        % Special solution for Wa = 0.9
        % a=(3^(1/2)*(3*10^(1/2))^(1/2))/9;
        % Wa=0.9;
        % b=(10 - 3*10^(1/2))^(1/2);
        % Wb=2/15;
        %
        % syms Wa real
        % syms Wb real
        % syms a real
        % syms b real
        % Wa = 0.95
        % ss=solve(8*Wa+6*Wb-8,8*a^2*Wa+2*b^2*Wb-8/3,8*a^4*Wa-8/9)
        % ss.a,  ss.b, ss.Wb
        % a=(57^(1/2)*(2*95^(1/2))^(1/2))/57;
        % Wa=0.95;
        % b=(20 - 2*95^(1/2))^(1/2);
        % Wb=1/15;
        %
        % syms Wb real
        % syms a real
        % syms b real
        % Wa = 0.3636
        % ss=solve(8*Wa+6*Wb-8,8*a^2*Wa+2*b^2*Wb-8/3,8*a^4*Wa+2*b^4*Wb-8/5)
        % ss.a, ss.Wa, ss.b, ss.Wb
        %
        % a=(303^(1/2)*(3408750/2159 - (15*1094305710^(1/2))/2159)^(1/2))/909;
        % Wa=0.3636;
        % b=((15*1094305710^(1/2))/3434969 + 1250/2159)^(1/2);
        % Wb=  1591/1875;
        %
        param_coords = [];
        weights = [];
    end
    
    methods
        
        function self = hex_14pt_rule (Parameters)
            self.param_coords = ...
                [self.a*[-1, -1, -1];...
                self.a*[+1, -1, -1];...
                self.a*[+1, +1, -1];...
                self.a*[-1, +1, -1];...
                self.a*[-1, -1, +1];...
                self.a*[+1, -1, +1];...
                self.a*[+1, +1, +1];...
                self.a*[-1, +1, +1];...
                self.b*[-1, 0, 0];...
                self.b*[+1, 0, 0];...
                self.b*[0, -1, 0];...
                self.b*[0, +1, 0];...
                self.b*[0, 0, -1];...
                self.b*[0, 0, +1];...
                ];
            self.weights =[self.Wa*ones(1,8),self.Wb*ones(1,6)];
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
