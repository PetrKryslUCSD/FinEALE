classdef  fe_set_T4 < fe_set_3_manifold
    % T4 (tetrahedron with four nodes)  finite element (FE) set class.
    %
    
    
    properties (Constant, GetAccess = public)
    end
    
    methods % constructor
        
        function self = fe_set_T4(Parameters)
            % Constructor.
            % Parameters:
            %            same as fe_set_3_manifold.
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin <1
                Parameters.nfens=4;
            end
            Parameters.nfens=4;
            self = self@fe_set_3_manifold(Parameters);
        end
        
    end
    
    
    methods % concrete methods
        
        function val = bfun(self,param_coords)
            % Evaluate the basis function matrix for four-node tetrahedron.
            %
            % function val = bfun(self,param_coords)
            %
            %   Call as:
            %      N = bfun (g, pc)
            %   where
            %      g=FE set,
            %      pc=parametric coordinates, -1 < pc(j) < 1, length(pc)=3.
            %
            val = [(1 - param_coords(1) - param_coords(2) - param_coords(3));...
                param_coords(1); ...
                param_coords(2); ...
                param_coords(3)];
        end
        
        function val = bfundpar (self, param_coords)
            % Evaluate the derivatives of the basis functions.
            %
            % function val = bfundpar (self, param_coords)
            %
            % Evaluate the derivatives of the basis function matrix for four-node
            % tetrahedron.
            % Returns an array of NFENS rows, and DIM columns, where
            %    NFENS=number of nodes, and
            %    DIM=number of spatial dimensions.
            % Call as:
            %    Nder = bfundpar(g, pc)
            % where g=FE set
            %       pc=parametric coordinates, -1 < pc(j) < 1, length(pc)=3.
            %
            val = [-1 -1 -1; ...
                +1  0  0; ...
                0 +1  0; ...
                0  0 +1];
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
            %    geom=geometry field
            % with optional fields
            % facecolor = color of the faces, solid
            % colorfield =field with vertex colors
            % only one of the facecolor and colorfield may be supplied
            %  shrink = shrink factor
            %
            % these are all defined in clockwise order, which is what Matlab expects
            context.faces=[1, 2, 3; 2, 4, 3; 3, 4, 1; 4, 2, 1];
            context.edges= [1,2;2,3;3,1;4,1;4,2;4,3];
            node_pc=[0,0,0;1,0,0;0,1,0;0,0,1];
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
            % Draw isosurface in a three-dimensional manifold finite element.
            %
            % function   draw_isosurface(self, gv, context)
            %
            % Input arguments
            % self = descendent of fe_set_3_manifold
            % gv = graphic viewer
            % context = struct
            % with mandatory fields
            %    x=reference geometry field
            %    u=displacement field
            % with optional fields
            %    scalarfield =field with scalar vertex data,
            %    isovalue = value of the isosurface
            %    color = what color should the isosurface be?  Default: [1,1,1]/2;
            %    edgealpha= transparency value
            %    linewidth=line width (default 2)
            conns = self.conn; % connectivity
            hconns=zeros(size(conns,1),8);
            hconns(:,1:4)=conns;
            hconns(:,5)=hconns(:,4);
            hconns(:,6)=hconns(:,4);
            hconns(:,7)=hconns(:,4);
            hconns(:,8)=hconns(:,4);
            hconns(:,4)=hconns(:,3);
            [f,v]= isosurf_faces(fe_set_H8, hconns, context.x, context.u, context.scalarfield, context.isovalue, 3);
            
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
        
        
        function val = in_parametric(self,param_coords)
            % Are given parametric coordinates inside the element parametric domain?
            %
            % function val = in_parametric(self,param_coords)
            %
            %  The parametric domain is the standard triangle.
            %
            % Returns a Boolean: is the point inside, true or false?
            %
            s =intersect(find(param_coords>=0),find(param_coords<=1));
            val = (length(s)==length(param_coords))&&(sum(param_coords)<=1);
        end
        
        function val =boundary_conn(self)
            % Get boundary connectivity.
            conn =self.conn;
            val = [conn(:,[1, 3, 2]);conn(:,[1, 2, 4]);conn(:,[2, 3, 4]);conn(:,[1, 4, 3])];
        end
        
        function val =boundary_fe(self)
            % Get the constructor of the class of the  boundary finite element.
            val = @fe_set_T3;
        end
        
    end
    
end


