classdef  fe_set_Q4 < fe_set_2_manifold
    % Q4 (quadrilateral with four nodes)  finite element (FE) set class.
    %
    %
    
    
    properties (Constant, GetAccess = public)
    end
    
    methods % constructor
        
        function self = fe_set_Q4(Parameters)
            % Constructor.
            % Parameters:
            % Parameters: same as fe_set_2_manifold.
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin <1
                Parameters.nfens=4;
            end
            Parameters.nfens=4;
            self = self@fe_set_2_manifold(Parameters);
        end
        
    end
    
    
    methods % concrete methods
        
        function val = bfun(self,param_coords)
            % Evaluate the basis function matrix for an 4-node quadrilateral.
            %
            % function val = bfun(self,param_coords)
            %
            %   Call as:
            %      N = bfun (g, pc)
            %   where
            %      g=FE set,
            %      pc=parametric coordinates, -1 < pc(j) < 1, length(pc)=3.
            %
            one_minus_xi = (1 - param_coords(1));
            one_plus_xi  = (1 + param_coords(1));
            one_minus_eta = (1 - param_coords(2));
            one_plus_eta  = (1 + param_coords(2));
            
            val = [0.25 * one_minus_xi * one_minus_eta;
                0.25 * one_plus_xi  * one_minus_eta;
                0.25 * one_plus_xi  * one_plus_eta;
                0.25 * one_minus_xi * one_plus_eta];
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
            xi  =  param_coords(1);
            eta =  param_coords(2);
            val = [ ...
                -(1. - eta) * 0.25, -(1. - xi) * 0.25;
                (1. - eta) * 0.25, -(1. + xi) * 0.25;
                (1. + eta) * 0.25, (1. + xi) * 0.25;
                -(1. + eta) * 0.25, (1. - xi) * 0.25];
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
            context.faces=[1 4 3 2];
            context.edges= [1,2,3,4,1];
            node_pc=[-1,-1;1,-1;1,1;-1,1];
            if isfield(context,'shrink')
                if (context.shrink~=1)
                    c=ones(size(node_pc,1),1)*(sum(node_pc)/size(node_pc,1));
                    pc=c+context.shrink*(node_pc-c);
                    for j=1:size(pc,1)
                        context.shrunk_pc_N{j}=bfun(self,pc(j,:));
                    end
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
            val = [conn(:,1:2);conn(:,2:3);conn(:,3:4);conn(:,[4,1])];
        end
        
        function val =boundary_fe(self)
            % Get the constructor of the class of the  boundary finite element.
            val = @fe_set_L2;
        end
        
    end
    
end


