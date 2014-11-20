classdef  fe_set_2_manifold < fe_set
    % Finite element (FE) set class.
    %
    % The set encloses some part of the computational domain.
    %
    % Class represents a two-manifold finite element set, which encloses some part
    % of the computational domain.
    % This class is the base class for surface-like finite elements, i.e.
    % every usable 2-D finite element has to be derived from it.
    %
    
    
    properties (Constant, GetAccess = public)
        dim=2;% manifold dimension (surface)
    end
    
    properties (GetAccess = public, SetAccess = protected)
        other_dimension = [];% the other dimension (in this case, cross-section area)
        axisymm = [];% is this an axially symmetric setting?
    end
    
    properties (Hidden, GetAccess = private, SetAccess = private)
        other_dimension_constant=[];
    end
    
    methods % constructor
        
        function self = fe_set_2_manifold(Parameters)
            % Constructor.
        % Parameters:
            %   same as fe_set.
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin <1
                Parameters = struct([]);
            end
            self = self@fe_set(Parameters);
            self.other_dimension = 1.0;
            if isfield(Parameters,'other_dimension')
                self.other_dimension = Parameters.other_dimension;
            end
            self.other_dimension_constant= false;
            if (strcmp(class(self.other_dimension),'double'))
                self.other_dimension_constant = true;
            end
            self.axisymm=false;
            if isfield(Parameters,'axisymm')
                self.axisymm = Parameters.axisymm;
            end
        end
        
    end
    
    
    methods % concrete methods
        
        function A= eval_other_dimension(self, conn, N, x)
            % Evaluate the other dimension (thickness).
            %
            % function A=eval_other_dimension(self, conn, N, x)
            %
            % Evaluate the other dimension (thickness) of the element at given
            % parametric coordinates, or at any given spatial coordinate. Arguments in
            % the same as for Jacobian.
            if (self.other_dimension_constant)
                A=self.other_dimension;
            else
                A=feval(self.other_dimension,conn, N, x);
            end
        end
        
        function Jac = Jacobian_surface(self, conn, N, J, x)
            % Evaluate the surface Jacobian.
            %
            % function Jac = Jacobian_surface(self, conn, N, J, x)
            %
            %   where
            %      self=FE set,
            %      conn=connectivity of a single cell in which the Jacobian is to be evaluated
            %      N=values of the basis functions, size(N) = [nbfuns,1]
            %      J = Jacobian matrix
            %      x=array of spatial coordinates of the nodes, size(x) = [nbfuns,dim]
            %
            [sdim, ntan] = size(J);
            if     ntan==2 % 2-D Set
                if sdim==ntan
                    Jac = J(1,1)*J(2,2)-J(2,1)*J(1,2);% isdet(J);% Compute the Jacobian
                else
                    Jac = norm(skewmat(J(:,1))*J(:,2));
                end
            else
                error('Got an incorrect size of tangents');
            end
        end
        
        
        function Jac = Jacobian_volume(self, conn, N, J, x)
            % Evaluate the volume Jacobian.
            %
            % function Jac = Jacobian_volume(self, conn, N, J, x)
            %
            %   For the two-dimensional cell, the volume Jacobian is
            %       (i) the product of the 2-D Jacobian and the other dimension (Units
            %       of length);
            %       or, when used as axially symmetric
            %       (ii) the product of the 2-D Jacobian and the circumference of
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
                Jac = Jacobian_surface(self, conn, N, J, x)*2*pi*xyz(1);
            else
                Jac = Jacobian_surface(self, conn, N, J, x)*eval_other_dimension(self, conn, N, x);
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
                otherwise
                    error('Wrong dimension');
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
            Jac = Jacobian_surface(self, conn, N, J, x);
        end
        
        function draw (self, gv, context)
            % Produce a graphic representation of the two-manifold FE set.
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
            ef=context.edges;
            f=context.faces;
            conns = self.conn; % connectivity
            if (isempty(conns))
                return; % nothing to be done
            end
            shrink  = 1.0;
            if isfield(context,'shrink')
                shrink =context.shrink;
            end
            if (shrink==1)
                nn =max(max(conns));
                x = gather_values(context.x, 1:nn); % coordinates of nodes
                u = gather_values(context.u, 1:nn); % coordinates of nodes
                xu=x+u;
                if (isfield(context,'colorfield'))
                    colors =gather_values(context.colorfield, 1:nn); % coordinates of nodes
                    context.colors=colors;
                    context.facecolor = 'none';
                end
                context1 = context;
                context1.edgecolor ='none';
                for i=1:size(f,1)
                  draw_polygon (gv, xu, conns(:,f(i,:)), context1);
                end
                if (~isfield(context,'edgecolor')) || (~strcmp(context.edgecolor,'none'))
                    draw_polyline (gv, xu, conns(:,ef), context);
                end
            else
                nfens=size(conns,2);
                for ac=1:size(conns,1)
                    conn =conns(ac,:);
                    x = gather_values(context.x, conn); % coordinates of nodes
                    u = gather_values(context.u, conn); % coordinates of nodes
                    xu=x+u;
                    xus=xu;
                    for j=1:length(context.shrunk_pc_N)
                        xus(j,:)=context.shrunk_pc_N{j}'*xu;
                    end
                    xu=xus;
                    if (isfield(context,'colorfield'))
                        colors =gather_values(context.colorfield, conn); % coordinates of nodes
                        context.colors=colors;
                        context.facecolor = 'none';
                    else
                        context.facecolor = 'red';
                        if (isfield(context,'facecolor'))
                            context.facecolor = context.facecolor;
                        end
                    end
                    if (isfield(context,'facealpha'))
                        context.facealpha = context.facealpha;
                    end
                    context1 = context;
                    context1.edgecolor ='none';
                    for i=1:size(f,1)
                        draw_polygon (gv, xu, f(i,:), context1);
                    end
                    if (~isfield(context,'edgecolor')) || (~strcmp(context.edgecolor,'none'))
                        draw_polyline (gv, xu, ef, context);
                    end
                end
            end
        end
        
        
        
    end
    
end


