classdef simpson_1_3_rule
% Class of Simpson 1/3 quadrature rule, on the interval -1 <=x<= +1.
%
% 


    properties
        dim = [];%  number of space dimensions (1,2,3)
        param_coords = [];% parametric coordinates of the quadrature points
        weights = [];  %  weights of the quadrature points
    end
    
    methods
    
        function self = simpson_1_3_rule (Parameters)
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
                    self.param_coords = [-1, 0, 1]';
                    self.weights = [1, 4, 1]/6*2;
                case 2
                    self.param_coords = [-1,-1; 1,-1; 1, 1; -1, 1; 0,-1; 1,0; 0,1; -1,0; 0,0];
                    w= 4*([1, 1, 4]/6)'*([1, 1, 4]/6);
                    self.weights = w([1     2     5     4     3     8     6     7     9]);
                case 3
                    p(1:4,1:3)= [-1,-1,-1;
                        1,-1,-1;
                        1,1,-1;
                        -1,1,-1];
                    p(5:8,1:3)= [-1,-1,1;
                        1,-1,1;
                        1,1,1;
                        -1,1,1];
                    p(9:12,1:3)= [0,-1,-1;
                        1,0,-1;
                        0,1,-1;
                        -1,0,-1];
                    p(13:16,1:3)= [0,-1,1;
                        1,0,1;
                        0,1,1;
                        -1,0,1];
                    p(17:20,1:3)= [-1,-1,0;
                        1,-1,0;
                        1,1,0;
                        -1,1,0];
                    p(21:26,1:3)= [0,0,-1;
                        0,-1,0;
                        1,0,0;
                        0,1,0;
                        -1,0,0;
                        0,0,1];
                    p(27,1:3)= [0,0,0];
                    self.param_coords = p;
                    w= 8*([1, 1, 4]/6)'*([1, 1, 4]/6);
                    W =zeros(3, 3, 3);
                    W(:,:,1)=w*(1/6);
                    W(:,:,2)=w*(1/6);
                    W(:,:,3)=w*(4/6);
                    self.weights = W([1     2     5     4    10    11    14    13     3     8     6     7    12    17    15    16    19    20    23   22     9    21    26    24    25    18    27]);
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


