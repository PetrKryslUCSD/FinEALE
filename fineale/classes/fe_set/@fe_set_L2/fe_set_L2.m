classdef  fe_set_L2 < fe_set_1_manifold
% L2 (curve segment with two nodes)  finite element (FE) set class.
%
%

    
    properties (Constant, GetAccess = public)
    end
    
    methods % constructor
        
        function self = fe_set_L2(Parameters)
            % Constructor.
            % Parameters:
            %           same as fe_set_1_manifold.
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin <1
                Parameters.nfens=2;
            end
            Parameters.nfens=2;
            self = self@fe_set_1_manifold(Parameters);
        end
        
    end
    
    
    methods % concrete methods
        
        function val = bfun(self,param_coords)
            % Evaluate the basis function matrix for an 2-node line element (bar).
            %
            % function val = bfun(self,param_coords)
            %
            %   Call as:
            %      N = bfun (g, pc)
            %   where
            %      g=FE set,
            %      pc=parametric coordinates, -1 < pc(j) < 1, length(pc)=3.
            %
            val = [(1 - param_coords(1)); (1 + param_coords(1))] / 2;
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
            Nder = [-1; +1]/2;
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
            context.edges= [1,2];
            node_pc=[-1,1];
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
            % edgecolor = 'Black';
            % if (isfield(context,'edgecolor'))
                % edgecolor =context.edgecolor;
            % end
            % if (isfield(context,'colorfield'))
            % else
                % facecolor = 'red';
                % if (isfield(context,'facecolor'))
                    % facecolor = context.facecolor;
                % end
                % edgecolor = 'Black';
                % if (isfield(context,'edgecolor'))
                    % edgecolor = context.edgecolor;
                % end
                % options=struct ('facecolor',facecolor,'edgecolor',edgecolor);
            % end
            % conns = self.conn; % connectivity
            % for ac=1:size(conns,1)
                % conn =conns(ac,:);
                % x = gather_values(context.x, conn); % coordinates of nodes
                % u = gather_values(context.u, conn); % coordinates of nodes
                % xu=x+u;
                % if isfield(context,'shrink')
                    % c=sum(xu,1)/2;
                    % xu=(1-context.shrink)*ones(2,1)*c+context.shrink*xu;
                % end
                % if (isfield(context,'colorfield'))
                    % colors =gather_values(context.colorfield, conn'); % coordinates of nodes
                    % options=struct ('colors',colors);
                % end
                % A1=other_dimension (self, -1);
                % A2=other_dimension (self, +1);
                % draw_cylinder(gv, xu(1,:),xu(2,:),sqrt(A1/pi),sqrt(A2/pi), context);
            % end
        end
        
        function val =boundary_conn(self)
            % Get boundary connectivity.
            val = [self.conn(:,1);self.conn(:,2)];
        end
        
        function val =boundary_fe(self)
            % Get the class of the  boundary finite element.
            val = @fe_set_P1;
        end
        
    end
    
end

