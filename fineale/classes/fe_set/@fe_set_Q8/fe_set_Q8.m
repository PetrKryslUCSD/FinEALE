classdef  fe_set_Q8 < fe_set_2_manifold
    % Q8 (quadrilateral with eight nodes)  finite element (FE) set class.
    %
    %
    
    
    properties (Constant, GetAccess = public)
    end
    
    methods % constructor
        
        function self = fe_set_Q8(Parameters)
            % Constructor.
            % Parameters:
            %            same as fe_set_2_manifold.
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin <1
                Parameters.nfens=8;
            end
            Parameters.nfens=8;
            self = self@fe_set_2_manifold(Parameters);
        end
        
    end
    
    
    methods % concrete methods
        
        function val = bfun(self,param_coords)
            % Evaluate the basis function matrix for an eight-node quadrilateral.
            %
            % function val = bfun(self,param_coords)
            %
            %   Call as:
            %      N = bfun (g, pc)
            %   where
            %      g=FE set,
            %      pc=parametric coordinates, -1 < pc(j) < 1, length(pc)=3.
            %
            xim    = (-1 + param_coords(1));
            etam   = (-1 + param_coords(2));
            xip    = (1 + param_coords(1));
            etap   = (1 + param_coords(2));
            val = [ -1.0/4*xim*etam*(1+param_coords(1)+param_coords(2));
                1.0/4*xip*etam*(1-param_coords(1)+param_coords(2));
                -1.0/4*xip*etap*(1-param_coords(1)-param_coords(2));
                1.0/4*xim*etap*(1+param_coords(1)-param_coords(2));
                1.0/2*xim*xip*etam;
                -1.0/2*etam*etap*xip;
                -1.0/2*xim*xip*etap;
                1.0/2*etam*etap*xim];
            return;
        end
        
        function val = bfundpar (self, param_coords)
            % Evaluate the derivatives of the basis function matrix.
            %
            % function val = bfundpar (self, param_coords)
            %
            % Returns an array of NFENS rows, and DIM columns, where
            %    NFENS=number of nodes, and
            %    DIM=number of spatial dimensions.
            % Call as:
            %    Nder = bfundpar(g, pc)
            % where g=FE set
            %       pc=parametric coordinates, -1 < pc(j) < 1, length(pc)=3.
            % Column j holds derivatives of the basis functions with respect
            % to the parametric coordinate j. Row i holds derivatives of the basis
            % function i with respect to all the parametric coordinates in turn.
            %
            xi =param_coords(1);
            eta =param_coords(2);
            n_xi= [1.0/4*(1-eta)*(1+xi+eta)-1.0/4*(1-xi)*(1-eta);
                -1.0/4*(1-eta)*(1-xi+eta)+1.0/4*(1+xi)*(1-eta);
                -1.0/4*(1+eta)*(1-xi-eta)+1.0/4*(1+xi)*(1+eta);
                1.0/4*(1+eta)*(1+xi-eta)-1.0/4*(1-xi)*(1+eta);
                -1.0/2*(1+xi)*(1-eta)+1.0/2*(1-xi)*(1-eta);
                1.0/2*(1-eta)*(1+eta);
                -1.0/2*(1+xi)*(1+eta)+1.0/2*(1-xi)*(1+eta);
                -1.0/2*(1-eta)*(1+eta)];
            n_eta = [1.0/4*(1-xi)*(1+xi+eta)-1.0/4*(1-xi)*(1-eta);
                1.0/4*(1+xi)*(1-xi+eta)-1.0/4*(1+xi)*(1-eta);
                -1.0/4*(1+xi)*(1-xi-eta)+1.0/4*(1+xi)*(1+eta);
                -1.0/4*(1-xi)*(1+xi-eta)+1.0/4*(1-xi)*(1+eta);
                -1.0/2*(1-xi)*(1+xi);
                -1.0/2*(1+xi)*(1+eta)+1.0/2*(1+xi)*(1-eta);
                1.0/2*(1-xi)*(1+xi);
                -1.0/2*(1-xi)*(1+eta)+1.0/2*(1-xi)*(1-eta)];
            val = [n_xi,n_eta];
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
            %    x=reference geometry field
            %    u=displacement field
            % with optional fields
            %    facecolor = color of the faces, solid
            %    colorfield =field with vertex colors
            %       only one of the facecolor and colorfield may be supplied
            %    shrink = shrink factor
            %
            context.faces=[1,5,8;2,6,5;3,7,6;4,8,7;5,6,7;7,8,5];
            context.edges= [1,5,2,6,3,7,4,8,1];
            node_pc=[-1,-1;1,-1;1,1;-1,1;0,-1;1,0;0,1;-1,0];
            if isfield(context,'shrink')
                c=ones(size(node_pc,1),1)*(sum(node_pc)/size(node_pc,1));
                pc=c+context.shrink*(node_pc-c);
                for j=1:size(pc,1)
                    context.shrunk_pc_N{j}=bfun(self,pc(j,:));
                end
            end
            draw@fe_set_2_manifold (self, gv, context);
        end
        
        function   draw_isosurface(self, gv, context)
            % Draw isosurface on a two-dimensional manifold FE set  (level curves, really).
            %
            % function   draw_isosurface(self, gv, context)
            %
            % Input arguments
            % self = descendent of fe_set_2_manifold
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
            conns = get(self, 'conn'); % connectivity
            [f,v]= isosurf_edges([conns(:, [1,2,3]);conns(:, [3,4,1])], context.x, context.u, context.scalarfield, context.isovalue);
            
            if (isfield(context,'color'))
                color =context.color;
            else
                color =[1,1,1]/2;
            end
            
            edgealpha = 1;
            if (isfield(context,'edgealpha'))
                edgealpha = context.edgealpha;
            end
            linewidth = 2;
            if (isfield(context,'linewidth'))
                linewidth = context.linewidth;
            end
            for j=1:size(f, 1 )
                patch('Faces',f(j,:),'Vertices',v,'FaceColor','none','EdgeColor',color,'linewidth',linewidth,'EdgeAlpha', edgealpha);
            end
        end
        
        function val =boundary_conn(self)
            % Get boundary connectivity.
            conn =self.conn;
            val = [conn(:,[1,2,5]);conn(:,[2,3,6]);conn(:,[3,4,7]);conn(:,[4,1,8])];
        end
        
        function val =boundary_fe(self)
            % Get the constructor of the class of the  boundary finite element.
            val = @fe_set_L3;
        end
        
    end
    
end


