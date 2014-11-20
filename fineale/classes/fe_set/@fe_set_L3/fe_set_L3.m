classdef  fe_set_L3 < fe_set_1_manifold
    % L3 (curve segment with three nodes)  finite element (FE) set class.
    %
    %
    
    
    
    properties (Constant, GetAccess = public)
    end
    
    methods % constructor
        
        function self = fe_set_L3(Parameters)
            % Constructor.
            % Parameters:
            % Parameters: same as fe_set_1_manifold.
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin <1
                Parameters.nfens=3;
            end
            Parameters.nfens=3;
            self = self@fe_set_1_manifold(Parameters);
        end
        
    end
    
    
    methods % concrete methods
        
        function val = bfun(self,param_coords)
            % Evaluate the basis function matrix for an 3-node line element (bar).
            %
            % function val = bfun(self,param_coords)
            %
            %   Call as:
            %      N = bfun (g, pc)
            %   where
            %      g=FE set,
            %      pc=parametric coordinates, -1 < pc(j) < 1, length(pc)=3.
            %
            xi=param_coords(1);
            val = [(xi-1)*xi/2;  (xi+1)*xi/2;  -(xi+1)*(xi-1)];
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
            %    Nder = bfundpar(g, pc)
            % where g=FE set
            %       pc=parametric coordinates, -1 < pc(j) < 1, length(pc)=dim.
            %
            xi=param_coords(1);
            Nder = [(xi-1/2); (xi+1/2); -2*xi];
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
            context.edges= [1,3,2];
            node_pc=[-1,0,1];
            if isfield(context,'shrink')
                if (context.shrink~=1)
                    c=ones(size(node_pc,1),1)*(sum(node_pc)/size(node_pc,1));
                    pc=c+context.shrink*(node_pc-c);
                    for j=1:size(pc,1)
                        context.shrunk_pc_N{j}=bfun(self,pc(j,:));
                    end
                end
            end
            draw@fe_set_1_manifold (self, gv, context);
            % conn = get(self, 'conn'); % connectivity
            % shrink = 1.0;
            % if isfield(context,'shrink')
                % shrink =context.shrink;
            % end
            % if (shrink~= 1.0)
                % xyz1 = map_to_xyz(self,-shrink,xu);
                % xyz2 = map_to_xyz(self,+shrink,xu);
                % xu(1,:)=xyz1;
                % xu(2,:)=xyz2;
            % end
            % if (isfield(context,'colorfield'))
                % colors =gather(context.colorfield, conn, 'values', 'noreshape'); % coordinates of nodes
                % options=struct ('colors',colors);
            % else
                % facecolor = 'red';
                % if (isfield(context,'facecolor'))
                    % facecolor = context.facecolor;
                % end
                % options=struct ('facecolor',facecolor);
            % end
            % A1=other_dimension (self, -shrink);
            % A2=other_dimension (self, +shrink);
            % A3=other_dimension (self, 0);
            % for j=1:size(conn,1)
                % x = gather(context.x, conn(j,:), 'values', 'noreshape'); % coordinates of nodes
                % u = gather(context.u, conn(j,:), 'values', 'noreshape'); % coordinates of nodes
                % xu=x+u;
                % draw_cylinder(gv, xu(1,:),xu(3,:),sqrt(A1/pi),sqrt(A3/pi), context);
                % draw_cylinder(gv, xu(3,:),xu(2,:),sqrt(A3/pi),sqrt(A2/pi), context);
            % end
        end
        
        
        function val =boundary_conn(self)
            % Get boundary connectivity.
            val = [self.conn(:,1);self.conn(:,2)];
        end
        
        function val =boundary_fe(self)
            % Get the constructor of the class of the  boundary finite element.
            val = @fe_set_P1;
        end
        
    end
    
end

