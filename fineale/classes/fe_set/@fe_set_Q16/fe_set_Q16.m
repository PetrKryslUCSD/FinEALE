classdef  fe_set_Q16 < fe_set_2_manifold
    % Q16 (quadrilateral with 16 nodes)  finite element (FE) set class.
    %
    %
    
    
    properties (Constant, GetAccess = public)
    end
    
    methods % constructor
        
        function self = fe_set_Q16(Parameters)
            % Constructor.
            % Parameters:
            %            same as fe_set_2_manifold.
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin <1
                Parameters.nfens=16;
            end
            Parameters.nfens=16;
            self = self@fe_set_2_manifold(Parameters);
        end
        
    end
    
    
    methods % concrete methods
        
        function N = bfun(self,param_coords)
            % Evaluate the basis function matrix for an 16-node quadrilateral.
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
            et=param_coords(2);
            N=[1/256*(xi-1)*(et-1)*(9*xi^2-1)*(9*et^2-1);
                -9/256*(3*xi-1)*(et-1)*(xi^2-1)*(9*et^2-1);
                9/256*(3*xi+1)*(et-1)*(xi^2-1)*(9*et^2-1);
                -1/256*(xi+1)*(et-1)*(9*xi^2-1)*(9*et^2-1);
                -9/256*(xi-1)*(3*et-1)*(9*xi^2-1)*(et^2-1);
                81/256*(3*xi-1)*(3*et-1)*(xi^2-1)*(et^2-1);
                -81/256*(3*xi+1)*(3*et-1)*(xi^2-1)*(et^2-1);
                9/256*(xi+1)*(3*et-1)*(9*xi^2-1)*(et^2-1);
                9/256*(xi-1)*(3*et+1)*(9*xi^2-1)*(et^2-1);
                -81/256*(3*xi-1)*(3*et+1)*(xi^2-1)*(et^2-1);
                81/256*(3*xi+1)*(3*et+1)*(xi^2-1)*(et^2-1);
                -9/256*(xi+1)*(3*et+1)*(9*xi^2-1)*(et^2-1);
                -1/256*(xi-1)*(et+1)*(9*xi^2-1)*(9*et^2-1);
                9/256*(3*xi-1)*(et+1)*(xi^2-1)*(9*et^2-1);
                -9/256*(3*xi+1)*(et+1)*(xi^2-1)*(9*et^2-1);
                1/256*(xi+1)*(et+1)*(9*xi^2-1)*(9*et^2-1)];
        end
        
        function N = bfundpar (self, param_coords)
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
            xi=param_coords(1);
            et=param_coords(2);
            N(1,:)=[1/256*(et-1)*(3*et-1)*(3*et+1)*(27*xi^2-18*xi-1),1/256*(-1-18*et+27*et^2)*(xi-1)*(3*xi-1)*(3*xi+1)];
            N(2,:)=[-9/256*(et-1)*(3*et-1)*(3*et+1)*(9*xi^2-3-2*xi),-9/256*(27*et^2-1-18*et)*(xi-1)*(3*xi-1)*(xi+1)];
            N(3,:)=[9/256*(et-1)*(3*et-1)*(3*et+1)*(-3+2*xi+9*xi^2),9/256*(27*et^2-1-18*et)*(xi-1)*(3*xi+1)*(xi+1)];
            N(4,:)=[-1/256*(et-1)*(3*et-1)*(3*et+1)*(27*xi^2-1+18*xi),-1/256*(27*et^2-1-18*et)*(3*xi-1)*(3*xi+1)*(xi+1)];
            N(5,:)=[-9/256*(et-1)*(3*et-1)*(et+1)*(27*xi^2-1-18*xi),-9/256*(9*et^2-3-2*et)*(xi-1)*(3*xi-1)*(3*xi+1)];
            N(6,:)=[81/256*(et-1)*(3*et-1)*(et+1)*(9*xi^2-3-2*xi),81/256*(-3-2*et+9*et^2)*(xi-1)*(3*xi-1)*(xi+1)];
            N(7,:)=[-81/256*(et-1)*(3*et-1)*(et+1)*(-3+2*xi+9*xi^2),-81/256*(-3-2*et+9*et^2)*(xi-1)*(3*xi+1)*(xi+1)];
            N(8,:)=[9/256*(et-1)*(3*et-1)*(et+1)*(27*xi^2+18*xi-1),9/256*(-3-2*et+9*et^2)*(3*xi-1)*(3*xi+1)*(xi+1)];
            N(9,:)=[9/256*(et-1)*(3*et+1)*(et+1)*(-1-18*xi+27*xi^2),9/256*(9*et^2-3+2*et)*(xi-1)*(3*xi-1)*(3*xi+1)];
            N(10,:)=[-81/256*(et-1)*(3*et+1)*(et+1)*(9*xi^2-2*xi-3),-81/256*(9*et^2-3+2*et)*(xi-1)*(3*xi-1)*(xi+1)];
            N(11,:)=[81/256*(et-1)*(3*et+1)*(et+1)*(2*xi+9*xi^2-3),81/256*(9*et^2+2*et-3)*(xi-1)*(3*xi+1)*(xi+1)];
            N(12,:)=[-9/256*(et-1)*(3*et+1)*(et+1)*(-1+27*xi^2+18*xi),-9/256*(9*et^2+2*et-3)*(3*xi-1)*(3*xi+1)*(xi+1)];
            N(13,:)=[-1/256*(3*et-1)*(3*et+1)*(et+1)*(-1+27*xi^2-18*xi),-1/256*(18*et+27*et^2-1)*(xi-1)*(3*xi-1)*(3*xi+1)];
            N(14,:)=[9/256*(3*et-1)*(3*et+1)*(et+1)*(9*xi^2-3-2*xi),9/256*(18*et+27*et^2-1)*(xi-1)*(3*xi-1)*(xi+1)];
            N(15,:)=[-9/256*(3*et-1)*(3*et+1)*(et+1)*(2*xi+9*xi^2-3),-9/256*(18*et-1+27*et^2)*(xi-1)*(3*xi+1)*(xi+1)];
            N(16,:)=[1/256*(3*et-1)*(3*et+1)*(et+1)*(27*xi^2-1+18*xi),1/256*(18*et-1+27*et^2)*(3*xi-1)*(3*xi+1)*(xi+1)];
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
            context.faces=[1,5,6,2;2,6,7,3;3,7,8,4;5,9,10,6;6,10,11,7;7,11,12,8;9,13,14,10;10,14,15,11;11,15,16,12];
            context.edges= [1,2,3,4,8,12,16,15,14,13,9,5,1];
            a=1/3;
            node_pc=[-1,-1;-a,-1;a,-1;1,-1;...
                -1,-a;-a,-a;a,-a;1,-a;...
                -1,+a;-a,+a;a,+a;1,+a;...
                -1,+1;-a,+1;a,+1;1,+1;...
                ];
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
            conns = self.conn; % connectivity
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
            val = [conn(:,[1,2,3,4]);conn(:,[4,8,12,16]);conn(:,[16,15,14,13]);conn(:,[13,9,5,1])];
        end
        
        function val =boundary_fe(self)
            % Get the constructor of the class of the  boundary finite element.
            val = @fe_set_L4;
        end
        
    end
    
end


