classdef gauss_rule
% Class of the Gauss rule.
%
% The rule is applicable for a tensor product of  intervals -1 <=x<= +1.
% 

    
    properties
        dim = [];% number of spatial dimensions (1,2,3)
        order = [];% order of the rule (number of quadrature points along the intervals -1 <=x<= +1)
        param_coords = [];% parametric coordinates of the quadrature points
        weights = [];  %  weights of the quadrature points
    end
    
    methods
        
        function self = gauss_rule (Parameters)
            % Constructor.
            % Parameters: 
            %   dim=number of  space dimensions
            %            3=3-D rule, 2=2-D  rule, 1=1-D rule (the default).
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin < 1
                return
            end
            self.dim = Parameters.dim;
            self.order = Parameters.order;
            self.param_coords=[];
            self.weights=[];
            switch self.order
                case 1
                    self.param_coords = [ 0 ];
                    self.weights  = [ 2 ];
                case 2
                    self.param_coords = [ -0.577350269189626 0.577350269189626 ];
                    self.weights  = [ 1 1 ];
                case 3
                    self.param_coords = [ -0.774596669241483  0  0.774596669241483 ];
                    self.weights  = [ 0.5555555555555556 0.8888888888888889 0.5555555555555556 ];
                case 4
                    self.param_coords = [ -0.86113631159405  -0.33998104358486   0.33998104358486   0.86113631159405];
                    self.weights  = [ 0.34785484513745   0.65214515486255   0.65214515486255   0.34785484513745];
                otherwise
                    [xx, ww] = gaussquad(self.order);
                    self.param_coords = xx';
                    self.weights  = ww';
            end
            switch self.dim
                case 1
                    self.param_coords = transpose (self.param_coords);
                    self.weights = self.weights;
                case 2
                    pc=self.param_coords;
                    w=self.weights;
                    self.param_coords=[];
                    self.weights=[];
                    for i=1:self.order
                        for j=1:self.order
                            self.param_coords = [ self.param_coords; [ pc(i) pc(j) ] ];
                            self.weights = [ self.weights; w(i)*w(j) ];
                        end
                    end
                case 3
                    pc=self.param_coords;
                    w=self.weights;
                    self.param_coords=[];
                    self.weights=[];
                    for i=1:self.order
                        for j=1:self.order
                            for k=1:self.order
                                self.param_coords = [ self.param_coords; [ pc(i) pc(j) pc(k) ] ];
                                self.weights = [ self.weights; w(i)*w(j)*w(k) ];
                            end
                        end
                    end
                otherwise
                    error('Unknown dimension!');
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
