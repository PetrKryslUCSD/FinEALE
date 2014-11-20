classdef  fe_set_H20 < fe_set_3_manifold
% H20 (hexahedron with 20 nodes)  finite element (FE) set class.
%
%

    
    properties (Constant, GetAccess = public)
    end
    
    methods % constructor
        
        function self = fe_set_H20(Parameters)
            % Constructor.
            % Parameters:
            %           same as fe_set_3_manifold.
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin <1
                Parameters.nfens=20;
            end
            Parameters.nfens=20;
            self = self@fe_set_3_manifold(Parameters);
        end
        
    end
    
    
    methods % concrete methods
        
        function val = bfun(self,param_coords)
            % Evaluate the basis function matrix for an 20-node brick.
            %
            % function val = bfun(self,param_coords)
            %
            %   Call as:
            %      N = bfun (g, pc)
            %   where
            %      g=gcell,
            %      pc=parametric coordinates, -1 < pc(j) < 1, length(pc)=3.
            %
            xim    = (-1 + param_coords(1));
            etam   = (-1 + param_coords(2));
            zetam  = (-1 + param_coords(3));
            xip    = (1 + param_coords(1));
            etap   = (1 + param_coords(2));
            zetap  = (1 + param_coords(3));
            val = [ 1.0/8*xim*etam*zetam*(2+param_coords(1)+param_coords(2)+param_coords(3));
                -1.0/8*xip*etam*zetam*(2-param_coords(1)+param_coords(2)+param_coords(3));
                1.0/8*xip*etap*zetam*(2-param_coords(1)-param_coords(2)+param_coords(3));
                -1.0/8*xim*etap*zetam*(2+param_coords(1)-param_coords(2)+param_coords(3));
                1.0/8*xim*etam*zetap*(-2-param_coords(1)-param_coords(2)+param_coords(3));
                -1.0/8*xip*etam*zetap*(-2+param_coords(1)-param_coords(2)+param_coords(3));
                1.0/8*xip*etap*zetap*(-2+param_coords(1)+param_coords(2)+param_coords(3));
                -1.0/8*xim*etap*zetap*(-2-param_coords(1)+param_coords(2)+param_coords(3));
                -1.0/4*xim*xip*etam*zetam;
                1.0/4*etam*etap*xip*zetam;
                1.0/4*xim*xip*etap*zetam;
                -1.0/4*etam*etap*xim*zetam;
                1.0/4*xim*xip*etam*zetap;
                -1.0/4*etam*etap*xip*zetap;
                -1.0/4*xim*xip*etap*zetap;
                1.0/4*etam*etap*xim*zetap;
                -1.0/4*zetam*zetap*xim*etam;
                1.0/4*zetam*zetap*xip*etam;
                -1.0/4*zetam*zetap*xip*etap;
                1.0/4*zetam*zetap*xim*etap];
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
            % where g=gcell
            %       pc=parametric coordinates, -1 < pc(j) < 1, length(pc)=3.
            %
            xim    = -(-1 + param_coords(1));
            etam   = -(-1 + param_coords(2));
            zetam  = -(-1 + param_coords(3));
            xip    = (1 + param_coords(1));
            etap   = (1 + param_coords(2));
            zetap  = (1 + param_coords(3));
            twoppp =(2+param_coords(1)+param_coords(2)+param_coords(3));
            twompp =(2-param_coords(1)+param_coords(2)+param_coords(3));
            twopmp =(2+param_coords(1)-param_coords(2)+param_coords(3));
            twoppm =(2+param_coords(1)+param_coords(2)-param_coords(3));
            twommp =(2-param_coords(1)-param_coords(2)+param_coords(3));
            twopmm =(2+param_coords(1)-param_coords(2)-param_coords(3));
            twompm =(2-param_coords(1)+param_coords(2)-param_coords(3));
            twommm =(2-param_coords(1)-param_coords(2)-param_coords(3));
            n_xi= [1.0/8*etam*zetam*twoppp-1.0/8*xim*etam*zetam;
                -1.0/8*etam*zetam*twompp+1.0/8*xip*etam*zetam;
                -1.0/8*etap*zetam*twommp+1.0/8*xip*etap*zetam;
                1.0/8*etap*zetam*twopmp-1.0/8*xim*etap*zetam;
                1.0/8*etam*zetap*twoppm-1.0/8*xim*etam*zetap;
                -1.0/8*etam*zetap*twompm+1.0/8*xip*etam*zetap;
                -1.0/8*etap*zetap*twommm+1.0/8*xip*etap*zetap;
                1.0/8*etap*zetap*twopmm-1.0/8*xim*etap*zetap;
                -1.0/4*xip*etam*zetam+1.0/4*xim*etam*zetam;
                1.0/4*etam*etap*zetam;
                -1.0/4*xip*etap*zetam+1.0/4*xim*etap*zetam;
                -1.0/4*etam*etap*zetam;
                -1.0/4*xip*etam*zetap+1.0/4*xim*etam*zetap;
                1.0/4*etam*etap*zetap;
                -1.0/4*xip*etap*zetap+1.0/4*xim*etap*zetap;
                -1.0/4*etam*etap*zetap;
                -1.0/4*zetam*zetap*etam;
                1.0/4*zetam*zetap*etam;
                1.0/4*zetam*zetap*etap;
                -1.0/4*zetam*zetap*etap];
            n_eta = [1.0/8*xim*zetam*twoppp-1.0/8*xim*etam*zetam;
                1.0/8*xip*zetam*twompp-1.0/8*xip*etam*zetam;
                -1.0/8*xip*zetam*twommp+1.0/8*xip*etap*zetam;
                -1.0/8*xim*zetam*twopmp+1.0/8*xim*etap*zetam;
                1.0/8*xim*zetap*twoppm-1.0/8*xim*etam*zetap;
                1.0/8*xip*zetap*twompm-1.0/8*xip*etam*zetap;
                -1.0/8*xip*zetap*twommm+1.0/8*xip*etap*zetap;
                -1.0/8*xim*zetap*twopmm+1.0/8*xim*etap*zetap;
                -1.0/4*xim*xip*zetam;
                -1.0/4*xip*etap*zetam+1.0/4*xip*etam*zetam;
                1.0/4*xim*xip*zetam;
                -1.0/4*xim*etap*zetam+1.0/4*xim*etam*zetam;
                -1.0/4*xim*xip*zetap;
                -1.0/4*xip*etap*zetap+1.0/4*xip*etam*zetap;
                1.0/4*xim*xip*zetap;
                -1.0/4*xim*etap*zetap+1.0/4*xim*etam*zetap;
                -1.0/4*zetam*zetap*xim;
                -1.0/4*zetam*zetap*xip;
                1.0/4*zetam*zetap*xip;
                1.0/4*zetam*zetap*xim];
            n_zeta= [1.0/8*xim*etam*twoppp-1.0/8*xim*etam*zetam;
                1.0/8*xip*etam*twompp-1.0/8*xip*etam*zetam;
                1.0/8*xip*etap*twommp-1.0/8*xip*etap*zetam;
                1.0/8*xim*etap*twopmp-1.0/8*xim*etap*zetam;
                -1.0/8*xim*etam*twoppm+1.0/8*xim*etam*zetap;
                -1.0/8*xip*etam*twompm+1.0/8*xip*etam*zetap;
                -1.0/8*xip*etap*twommm+1.0/8*xip*etap*zetap;
                -1.0/8*xim*etap*twopmm+1.0/8*xim*etap*zetap;
                -1.0/4*xim*xip*etam;
                -1.0/4*etam*etap*xip;
                -1.0/4*xim*xip*etap;
                -1.0/4*etam*etap*xim;
                1.0/4*xim*xip*etam;
                1.0/4*etam*etap*xip;
                1.0/4*xim*xip*etap;
                1.0/4*etam*etap*xim;
                -1.0/4*xim*etam*zetap+1.0/4*xim*etam*zetam;
                -1.0/4*xip*etam*zetap+1.0/4*xip*etam*zetam;
                -1.0/4*xip*etap*zetap+1.0/4*xip*etap*zetam;
                -1.0/4*xim*etap*zetap+1.0/4*xim*etap*zetam];
            Nder = [n_xi, n_eta, n_zeta];
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
            context.faces=[1,9,12,12;9,2,10,10;10,3,11,11;11,4,12,12;9,10,11,12;
                5,16,13,13;16,8,15,15;15,7,14,14;14,6,13,13;16,15,14,13;
                1,17,9,9;17,5,13,13;13,6,18,18;18,2,9,9;17,13,18,9;
                2,18,10,10;18,6,14,14;14,7,19,19;19,3,10,10;18,14,19,10;
                3,19,4,4;19,7,15,15;15,8,20,20;20,4,11,11;19,15,20,11;
                4,20,12,12;20,8,16,16;16,5,17,17;17,1,12,12;20,16,17,12];
            context.edges= [1,9,2,18,6,13,5;
                2,10,3,19,7,14,6;
                3,11,4,20,8,15,7;
                4,12,1,17,5,16,8];
            node_pc=[-1,-1,-1;+1,-1,-1;+1,+1,-1;-1,+1,-1;...
                -1,-1,+1;+1,-1,+1;+1,+1,+1;-1,+1,+1;...
                0,-1,-1;1,0,-1;0,1,-1;-1,0,-1;...
                0,-1,1;1,0,1;0,1,1;-1,0,1;...
                -1,-1,0;+1,-1,0;+1,+1,0;-1,+1,0];
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
            % Note: the nodes need to be specified counterclockwise on each
            % face when viewed from outside of the element
            val = [conn(:,[1, 4, 3, 2, 12, 11, 10, 9]);
                conn(:,[1, 2, 6, 5, 9, 18, 13, 17]);
                conn(:,[2, 3, 7, 6, 10, 19, 14, 18]);
                conn(:,[3, 4, 8, 7,  11, 20, 15, 19]);
                conn(:,[4, 1, 5, 8, 12, 17, 16, 20]);
                conn(:,[6, 7, 8, 5, 14, 15, 16, 13])];
        end
        
        function val =boundary_fe(self)
            % Get the constructor of the class of the  boundary finite element.
            val = @fe_set_Q8;
        end
        
    end
    
end


