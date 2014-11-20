classdef  fe_set_0_manifold < fe_set
    % Finite element (FE) set class of zero manifold dimension (point).
    %
    % Class represents a zero-manifold (point-like) finite element, which encloses some part
    % of the computational domain.
    % This class is the base class for point-like finite elements, i.e. every usable 0-D FE has to be
    % derived from it.
    %
    
    
    properties (Constant, GetAccess = public)
        dim=0;% manifold dimension (point)
    end
    
    properties (GetAccess = public, SetAccess = protected)
        other_dimension = [];% the other dimension of the manifold (in this case, volume)
        axisymm = [];% is the setting axially symmetric?
    end
    
    properties (Hidden, GetAccess = private, SetAccess = private)
        other_dimension_constant=[];
    end
    
    methods % constructor
        
        function self = fe_set_0_manifold(Parameters)
            % Constructor.
            % Parameters:
            %   same as fe_set.
            if nargin <1
                Parameters = struct([]);
            end
            self = self@fe_set(Parameters);
            self.other_dimension = 1.0;
            self.other_dimension_constant = true;
            if isfield(Parameters,'other_dimension')
                self.other_dimension_constant= false;
                if (strcmp(class(Parameters.other_dimension),'double'))
                    self.other_dimension = Parameters.other_dimension;
                    self.other_dimension_constant = true;
                else
                    self.other_dimension =Parameters.other_dimension;
                end
            end
            self.axisymm=false;
            if isfield(Parameters,'axisymm')
                self.axisymm = Parameters.axisymm;
            end
        end
        
    end
    
    
    methods % concrete methods
        
        function A= eval_other_dimension(self, conn, N, x)
            % Evaluate the other dimension (volume, thickness, length).
            %
            % function A=eval_other_dimension(self, conn, N, x)
            %
            % Evaluate the other dimension (volume, thickness, length) of the element at given
            % parametric coordinates, or at any given spatial coordinate. Arguments in
            % the same as for Jacobian.
            if (self.other_dimension_constant)
                A=self.other_dimension;
            else
                A=feval(self.other_dimension,conn, N, x);
            end
        end
        
        function Jac = Jacobian_curve(self, conn, N, J, x)
            % Evaluate the curve Jacobian.
            %
            % function Jac = Jacobian_curve(self, conn, N, J, x)
            %
            %
            %   For the zero-dimensional cell, the curve Jacobian is
            %       (i) the product of the 0-D Jacobian and the other dimension (units
            %       of length);
            %       or, when used as axially symmetric
            %       (ii) the product of the 0-D Jacobian and the circumference of
            %       the circle through the point pc
            %
            %   where
            %      self=FE set,
            %      conn=connectivity of a single cell in which the Jacobian is to be evaluated
            %      N=values of the basis functions, size(N) = [nbfuns,1]
            %      J = Jacobian matrix
            %      x=array of spatial coordinates of the nodes, size(x) = [nbfuns,dim]
            %
            if self.axisymm
                xyz =N'*x;
                Jac = 1*2*pi*xyz(1);
            else
                Jac = 1*eval_other_dimension(self, conn, N, x);
            end
        end
        
        function Jac = Jacobian_surface(self, conn, N, J, x)
            % Evaluate the surface Jacobian.
            %
            % function Jac = Jacobian_surface(self, conn, N, J, x)
            %
            %
            %   For the zero-dimensional cell, the surface Jacobian is
            %       (i) the product of the 0-D Jacobian and the other dimension (Units
            %       of length squared);
            %       or, when used as axially symmetric
            %       (ii) the product of the 0-D Jacobian and the circumference of
            %       the circle through the point pc and the other dimension (units of
            %       length)
            %
            %   where
            %      self=FE set,
            %      conn=connectivity of a single cell in which the Jacobian is to be evaluated
            %      N=values of the basis functions, size(N) = [nbfuns,1]
            %      J = Jacobian matrix
            %      x=array of spatial coordinates of the nodes, size(x) = [nbfuns,dim]
            %
            if self.axisymm
                xyz =N'*x;
                Jac = 1*2*pi*xyz(1)*eval_other_dimension(self, conn, N, x);
            else
                Jac = 1*eval_other_dimension(self, conn, N, x);
            end
        end
        
        
        function Jac = Jacobian_volume(self, conn, N, J, x)
            % Evaluate the volume Jacobian.
            %
            % function Jac = Jacobian_volume(self, conn, N, J, x)
            %
            %   For the zero-dimensional cell, the volume Jacobian is the product of the
            %    0-D Jacobian and the other dimension (units of length cubed).
            %
            %   where
            %      self=FE set,
            %      conn=connectivity of a single cell in which the Jacobian is to be evaluated
            %      N=values of the basis functions, size(N) = [nbfuns,1]
            %      J = Jacobian matrix
            %      x=array of spatial coordinates of the nodes, size(x) = [nbfuns,dim]
            %
            if self.axisymm
                xyz =N'*x;
                Jac = 1*2*pi*xyz(1)*eval_other_dimension(self, conn, N, x);
            else
                Jac = 1*eval_other_dimension(self, conn, N, x);
            end
        end
        
        function Jac = Jacobian_mdim(self, conn, N, J, x, m)
            % Evaluate the manifold Jacobian.
            %
            % function Jac = Jacobian_mdim(self, conn, N, J, x, m)
            %
            %
            %   For and m-dimensional cell, the manifold Jacobian is
            %       m=0: +1
            %       m=1: Jacobian_curve
            %       m=2: Jacobian_surface
            %       m=3: Jacobian_volume
            %
            %   where
            %      self=FE set,
            %      conn=connectivity of a single cell in which the Jacobian is to be evaluated
            %      N=values of the basis functions, size(N) = [nbfuns,1]
            %      J = Jacobian matrix
            %      x=array of spatial coordinates of the nodes, size(x) = [nbfuns,dim]
            %
            switch (m)
                case 3
                    Jac = Jacobian_volume(self, conn, N, J, x);
                case 2
                    Jac = Jacobian_surface(self, conn, N, J, x);
                case 1
                    Jac = Jacobian_curve(self, conn, N, J, x);
                case 0
                    Jac = 1;
                otherwise
                    Jac = Jacobian(self, conn, N, J, x);
            end
        end
        
        function Jac = Jacobian(self, conn, N, J, x)
            % Evaluate the manifold Jacobian.
            %
            % function Jac = Jacobian(self, conn, N, J, x)
            %
            %   For and m-dimensional cell, the manifold Jacobian is
            %       m=0: +1
            %       m=1: Jacobian_curve
            %       m=2: Jacobian_surface
            %       m=3: Jacobian_volume
            %
            %   where
            %      self=FE set,
            %      conn=connectivity of a single cell in which the Jacobian is to be evaluated
            %      N=values of the basis functions, size(N) = [nbfuns,1]
            %      J = Jacobian matrix
            %      x=array of spatial coordinates of the nodes, size(x) = [nbfuns,dim]
            %
            Jac = Jacobian_curve(self, conn, N, J, x);
        end
        
        function draw (self, gv, context)
            % Produce a graphic representation of the finite element.
            %
            % function draw (self, gv, context)
            %
            % Input arguments
            % self = self
            % gv = graphic viewer
            % context = struct
            % with mandatory fields
            %    x=reference geometry field
            %    u=displacement field
            % with optional fields
            %    facecolor = color of the faces, solid
            %    colorfield =field with vertex colors
            %       only one of the facecolor and colorfield may be supplied
            %    shrink = shrink factor
            %
            
            conns = self.conn; % connectivity
            if (isempty(conns))
                return; % nothing to be done
            end
            
            nn =max(max(conns));
            x = gather_values(context.x, conns); % coordinates of nodes
            u = gather_values(context.u, conns); % coordinates of nodes
            xu=x+u;
            if (isfield(context,'colorfield'))
                colors =gather_values(context.colorfield, 1:nn); % coordinates of nodes
                context.colors=colors;
                context.facecolor = 'none';
            end
            if isfield(context,'edgecolor') && (~strcmp(context.edgecolor,'none'))
                draw_marker(gv, xu, context);
            end
            
        end
        
        
        
    end
    
end


