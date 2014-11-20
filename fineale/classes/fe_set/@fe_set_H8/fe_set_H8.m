classdef  fe_set_H8 < fe_set_3_manifold
    % H8 (hexahedron with eight nodes)  finite element (FE) set class.
    %
    %
    
    
    properties (Constant, GetAccess = public)
    end
    
    methods % constructor
        
        function self = fe_set_H8(Parameters)
            % Constructor.
            % Parameters:
            %           same as fe_set_3_manifold.
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin <1
                Parameters.nfens=8;
            end
            Parameters.nfens=8;
            self = self@fe_set_3_manifold(Parameters);
        end
        
    end
    
    
    methods % concrete methods
        
        function val = bfun(self,param_coords)
            % Evaluate the basis function matrix for an 8-node brick.
            %
            % function val = bfun(self,param_coords)
            %
            %   Call as:
            %      N = bfun (g, pc)
            %   where
            %      g=FE set,
            %      pc=parametric coordinates, -1 < pc(j) < 1, length(pc)=3.
            %
            one_minus_xi    = (1 - param_coords(1));
            one_minus_eta   = (1 - param_coords(2));
            one_minus_theta = (1 - param_coords(3));
            one_plus_xi     = (1 + param_coords(1));
            one_plus_eta    = (1 + param_coords(2));
            one_plus_theta  = (1 + param_coords(3));
            val = [one_minus_xi*one_minus_eta*one_minus_theta;...
                one_plus_xi*one_minus_eta*one_minus_theta;...
                one_plus_xi*one_plus_eta*one_minus_theta;...
                one_minus_xi*one_plus_eta*one_minus_theta;...
                one_minus_xi*one_minus_eta*one_plus_theta;...
                one_plus_xi*one_minus_eta*one_plus_theta;...
                one_plus_xi*one_plus_eta*one_plus_theta;...
                one_minus_xi*one_plus_eta*one_plus_theta] / 8;
            return;
        end
        
        function Nder = bfundpar (self, param_coords)
            % Evaluate the derivatives of the basis function matrix.
            %
            % function Nder = bfundpar (self, param_coords)
            %
            % Returns an array of NFENS rows, and DIM columns, where
            %    NFENS=number of nodes, and
            %    DIM=number of spatial dimensions.
            % Call as:
            %    Nder = bfundpar(g, pc), the matrix of basis function gradients
            %      with respect to the parametric coordinates
            % where g=FE set
            %       pc=parametric coordinates, -1 < pc(j) < 1, length(pc)=3.
            %
            one_minus_xi    = (1 - param_coords(1));
            one_minus_eta   = (1 - param_coords(2));
            one_minus_theta = (1 - param_coords(3));
            one_plus_xi     = (1 + param_coords(1));
            one_plus_eta    = (1 + param_coords(2));
            one_plus_theta  = (1 + param_coords(3));
            Nder = [-one_minus_eta*one_minus_theta...
                one_minus_eta*one_minus_theta...
                one_plus_eta*one_minus_theta...
                -one_plus_eta*one_minus_theta...
                -one_minus_eta*one_plus_theta...
                one_minus_eta*one_plus_theta...
                one_plus_eta*one_plus_theta...
                -one_plus_eta*one_plus_theta;...
                -one_minus_xi*one_minus_theta...
                -one_plus_xi*one_minus_theta...
                one_plus_xi*one_minus_theta...
                one_minus_xi*one_minus_theta...
                -one_minus_xi*one_plus_theta...
                -one_plus_xi*one_plus_theta...
                one_plus_xi*one_plus_theta...
                one_minus_xi*one_plus_theta;...
                -one_minus_xi*one_minus_eta...
                -one_plus_xi*one_minus_eta...
                -one_plus_xi*one_plus_eta...
                -one_minus_xi*one_plus_eta...
                one_minus_xi*one_minus_eta...
                one_plus_xi*one_minus_eta...
                one_plus_xi*one_plus_eta...
                one_minus_xi*one_plus_eta]'/8;
            return;
        end
        
        function draw (self, gv, context)
            % Produce a graphic representation.
            %
            % function draw (self, gv, context)
            %
            %
            % Input arguments
            % self = self
            % gv = graphic viewer
            % context = struct
            % with mandatory fields
            %    x=reference geometry nodal_field
            %    u=displacement nodal_field
            % with optional fields
            %    facecolor = color of the faces, solid
            %    colorfield =nodal_field with vertex colors
            %       only one of the facecolor and colorfield may be supplied
            %    shrink = shrink factor
            %
            % these are all defined in clockwise order, which is what Matlab expects
            context.faces=[ 4     1     2     3
                2     1     5     6
                3     2     6     7
                4     3     7     8
                1     4     8     5
                6     5     8     7];
            context.edges= [1,2;2,3;3,4;4,1;5,6;6,7;7,8;8,5;1,5;2,6;3,7;4,8];
            node_pc=[-1,-1,-1;+1,-1,-1;+1,+1,-1;-1,+1,-1;...
                -1,-1,+1;+1,-1,+1;+1,+1,+1;-1,+1,+1;];
            if isfield(context,'shrink')
                c=ones(size(node_pc,1),1)*(sum(node_pc)/size(node_pc,1));
                pc=c+context.shrink*(node_pc-c);
                for j=1:size(pc,1)
                    context.shrunk_pc_N{j}=bfun(self,pc(j,:));
                end
            end
            draw@fe_set_3_manifold (self, gv, context);
        end
        
        function   draw_isosurface(self, gv, context)
            % Draw isosurface within a three-dimensional manifold finite element.
            %
            % function   draw_isosurface(self, gv, context)
            %
            % Input arguments
            % self = descendent of fe_set_3_manifold
            % gv = graphic viewer
            % context = struct
            % with mandatory fields
            %    x=reference geometry nodal_field
            %    u=displacement field
            % with optional fields
            %    scalarfield =field with scalar vertex data,
            %    isovalue = value of the isosurface
            %    color = what color should the isosurface be?  Default: [1,1,1]/2;
            %    edgealpha= transparency value
            %    linewidth=line width (default 2)
            conns = self.conn; % connectivity
            [f,v]= isosurf_faces(self, conns, context.x, context.u, context.scalarfield, context.isovalue, 3);
            
            if (isfield(context,'color'))
                color =context.color;
            else
                color =[1,1,1]/2;
            end
            
            facealpha = 1;
            if (isfield(context,'facealpha'))
                facealpha = context.facealpha;
            end
            patch('Faces',f,'Vertices',v,'FaceColor',color,'EdgeColor','none','FaceAlpha', facealpha);
        end
        
        
        
        function val =boundary_conn(self)
            % Get boundary connectivity.
            conn =self.conn;
            val = [conn(:,[1, 4, 3, 2]);
                conn(:,[1, 2, 6, 5]);
                conn(:,[2, 3, 7, 6]);
                conn(:,[3, 4, 8, 7]);
                conn(:,[4, 1, 5, 8]);
                conn(:,[6, 7, 8, 5])];
        end
        
        function val =boundary_fe(self)
            % Get the constructor of the class of the  boundary finite element.
            val = @fe_set_Q4;
        end
        
    end
    
end


