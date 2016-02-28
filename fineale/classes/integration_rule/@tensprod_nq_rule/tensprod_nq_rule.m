% Constructor of tensor-product nodal-quadrature rule.
% 
% function retobj = simplex_nq_rule (Parameters)
% 
%   Parameters=structure with fields
%      dim=number of  space dimensions
%            3=3-D rule (four-node tetrahedron), 2=2-D  rule (three-node triangle).
% 
%
classdef tensprod_nq_rule
    
    properties
        param_coords = [];
        weights = [];
    end
    
    methods
        
        function self = tensprod_nq_rule (Parameters)
            if nargin < 1
                return
            end
            dim = Parameters.dim;
            order =Parameters.order;
            if (order== 1)
                switch dim
                    case 1
                        self.param_coords = [-1; 1];
                        self.weights = [1, 1];
                    case 2
                        self.param_coords = [-1,-1; 1,-1; 1, 1; -1, 1];
                        self.weights = [1, 1, 1, 1];
                    case 3
                        self.param_coords = [-1,-1,-1; 1,-1,-1; 1, 1,-1;-1, 1,-1;-1,-1, 1; 1,-1, 1; 1, 1, 1;-1, 1, 1;];
                        self.weights = [1, 1, 1, 1, 1, 1, 1, 1];
                    otherwise
                        error('Unsupported Dimension');
                end
            elseif (order== 2)
                switch dim
                    case 1
                        self.param_coords = [-1; 1; 0];
                        self.weights = 2*[1, 1, 4]/6;
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
            elseif (self.order== 3)
                switch self.dim
                    case 1
                        self.param_coords = [-1; -1/3; 1/3; 1];
                        self.weights = 2/3*[1, 3, 3, 1]*3/8;
                    case 2
                        self.param_coords = [ -1, -1;
                            -0.33333, -1;
                            0.33333, -1;
                            1, -1;
                            -1, -0.33333;
                            -0.33333, -0.33333;
                            0.33333, -0.33333;
                            1, -0.33333;
                            -1, 0.33333;
                            -0.33333, 0.33333;
                            0.33333, 0.33333;
                            1, 0.33333;
                            -1, 1;
                            -0.33333, 1;
                            0.33333, 1;
                            1, 1;];
                        w= (2/3*[1, 3, 3, 1]*3/8)'*(2/3*[1, 3, 3, 1]*3/8);
                        self.weights = w((1:16));
                    case 3
                        p= [-1, -1, -1;
                            -0.33333, -1, -1;
                            0.33333, -1, -1;
                            1, -1, -1;
                            -1, -0.33333, -1;
                            -0.33333, -0.33333, -1;
                            0.33333, -0.33333, -1;
                            1, -0.33333, -1;
                            -1, 0.33333, -1;
                            -0.33333, 0.33333, -1;
                            0.33333, 0.33333, -1;
                            1, 0.33333, -1;
                            -1, 1, -1;
                            -0.33333, 1, -1;
                            0.33333, 1, -1;
                            1, 1, -1;
                            -1, -1, -0.33333;
                            -0.33333, -1, -0.33333;
                            0.33333, -1, -0.33333;
                            1, -1, -0.33333;
                            -1, -0.33333, -0.33333;
                            -0.33333, -0.33333, -0.33333;
                            0.33333, -0.33333, -0.33333;
                            1, -0.33333, -0.33333;
                            -1, 0.33333, -0.33333;
                            -0.33333, 0.33333, -0.33333;
                            0.33333, 0.33333, -0.33333;
                            1, 0.33333, -0.33333;
                            -1, 1, -0.33333;
                            -0.33333, 1, -0.33333;
                            0.33333, 1, -0.33333;
                            1, 1, -0.33333;
                            -1, -1, 0.33333;
                            -0.33333, -1, 0.33333;
                            0.33333, -1, 0.33333;
                            1, -1, 0.33333;
                            -1, -0.33333, 0.33333;
                            -0.33333, -0.33333, 0.33333;
                            0.33333, -0.33333, 0.33333;
                            1, -0.33333, 0.33333;
                            -1, 0.33333, 0.33333;
                            -0.33333, 0.33333, 0.33333;
                            0.33333, 0.33333, 0.33333;
                            1, 0.33333, 0.33333;
                            -1, 1, 0.33333;
                            -0.33333, 1, 0.33333;
                            0.33333, 1, 0.33333;
                            1, 1, 0.33333;
                            -1, -1, 1;
                            -0.33333, -1, 1;
                            0.33333, -1, 1;
                            1, -1, 1;
                            -1, -0.33333, 1;
                            -0.33333, -0.33333, 1;
                            0.33333, -0.33333, 1;
                            1, -0.33333, 1;
                            -1, 0.33333, 1;
                            -0.33333, 0.33333, 1;
                            0.33333, 0.33333, 1;
                            1, 0.33333, 1;
                            -1, 1, 1;
                            -0.33333, 1, 1;
                            0.33333, 1, 1;
                            1, 1, 1];
                        self.param_coords = p;
                        w= (2/3*[1, 3, 3, 1]*3/8)'*(2/3*[1, 3, 3, 1]*3/8);
                        W =zeros(4,4,4);
                        W(:,:,1)=w*(2/3*3/8);
                        W(:,:,2)=w*(2/3*3*3/8);
                        W(:,:,3)=w*(2/3*3*3/8);
                        W(:,:,4)=w*(2/3*3/8);
                        self.weights = W((1:64));
                    otherwise
                        error('Unsupported Dimension');
                end
            else
                error('Unsupported  order');
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
% 