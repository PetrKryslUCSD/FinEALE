classdef  fe_set_H27 < fe_set_3_manifold
% H27 (hexahedron with 27 nodes)  finite element (FE) set class.
%
%

    
    properties (Constant, GetAccess = public)
    end
    
    methods % constructor
        
        function self = fe_set_H27(Parameters)
            % Constructor.
            % Parameters:
            %           same as fe_set_3_manifold.
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin <1
                Parameters.nfens=27;
            end
            Parameters.nfens=27;
            self = self@fe_set_3_manifold(Parameters);
        end
        
    end
    
    
    methods % concrete methods
        
        function val = bfun(self,param_coords)
            % Evaluate the basis function matrix for an 27-node brick.
            %
            % function val = bfun(self,param_coords)
            %
            %   Call as:
            %      N = bfun (g, pc)
            %   where
            %      g=gcell,
            %      pc=parametric coordinates, -1 < pc(j) < 1, length(pc)=3.
            %
            %     p(1:4,1:3)= [-1,-1,-1;
            %         1,-1,-1;
            %         1,1,-1;
            %         -1,1,-1];
            %     p(5:8,1:3)= [-1,-1,1;
            %         1,-1,1;
            %         1,1,1;
            %         -1,1,1];
            %     p(9:12,1:3)= [0,-1,-1;
            %         1,0,-1;
            %         0,1,-1;
            %         -1,0,-1];
            %     p(13:16,1:3)= [0,-1,1;
            %         1,0,1;
            %         0,1,1;
            %         -1,0,1];
            %     p(17:20,1:3)= [-1,-1,0;
            %         1,-1,0;
            %         1,1,0;
            %         -1,1,0];
            %     p(21:26,1:3)= [0,0,-1;
            %         0,-1,0;
            %         1,0,0;
            %         0,1,0;
            %         -1,0,0;
            %         0,0,1];
            %     p(27,1:3)= [0,0,0];
            %     idx=[1     2     5     4    10    11    14    13     3     8     6     7    12    17    15    16    19    20    23   22     9    21    26    24    25    18    27];
            %     %         idx=[];
            %     %     for i= 1 : 27
            % %     param_coords=p(i,:);
            %     xi=param_coords(1);
            %     xis = [(xi-1)*xi/2;  (xi+1)*xi/2;  -(xi+1)*(xi-1)];
            %     eta=param_coords(2);
            %     etas = [(eta-1)*eta/2;  (eta+1)*eta/2;  -(eta+1)*(eta-1)];
            %     zeta=param_coords(3);
            %     zetas = [(zeta-1)*zeta/2;  (zeta+1)*zeta/2;  -(zeta+1)*(zeta-1)];
            %     N =zeros(3, 3, 3);
            %     xisetas=(xis*etas');
            %     N(:,:,1)=xisetas*zetas(1);
            %     N(:,:,2)=xisetas*zetas(2);
            %     N(:,:,3)=xisetas*zetas(3);
            %     %     N
            %     %         idx=[idx,find(N==1)]
            %     %     end
            %     val = N(idx)';
            xi=param_coords(1);
            eta=param_coords(2);
            zet=param_coords(3);
            x1 =(xi-1);
            y1=(eta-1);
            z1 =(zet-1);
            x2 =(xi+1);
            y2=(eta+1);
            z2 =(zet+1);
            val =[1/8*z1*zet*x1*xi*y1*eta
                1/8*z1*zet*x2*xi*y1*eta
                1/8*z1*zet*x2*xi*y2*eta
                1/8*z1*zet*x1*xi*y2*eta
                1/8*z2*zet*x1*xi*y1*eta
                1/8*z2*zet*x2*xi*y1*eta
                1/8*z2*zet*x2*xi*y2*eta
                1/8*z2*zet*x1*xi*y2*eta
                1/4*z1*zet*(-x2)*x1*y1*eta
                1/4*z1*zet*x2*xi*(-y2)*y1
                1/4*z1*zet*(-x2)*x1*y2*eta
                1/4*z1*zet*x1*xi*(-y2)*y1
                1/4*z2*zet*(-x2)*x1*y1*eta
                1/4*z2*zet*x2*xi*(-y2)*y1
                1/4*z2*zet*(-x2)*x1*y2*eta
                1/4*z2*zet*x1*xi*(-y2)*y1
                1/4*(-z2)*z1*x1*xi*y1*eta
                1/4*(-z2)*z1*x2*xi*y1*eta
                1/4*(-z2)*z1*x2*xi*y2*eta
                1/4*(-z2)*z1*x1*xi*y2*eta
                1/2*z1*zet*(-x2)*x1*(-y2)*y1
                1/2*(-z2)*z1*(-x2)*x1*y1*eta
                1/2*(-z2)*z1*x2*xi*(-y2)*y1
                1/2*(-z2)*z1*(-x2)*x1*y2*eta
                1/2*(-z2)*z1*x1*xi*(-y2)*y1
                1/2*z2*zet*(-x2)*x1*(-y2)*y1
                (-z2)*z1*(-x2)*x1*(-y2)*y1];
            return;
        end
        
        
        function val = bfundpar (self, param_coords)
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
            %     idx=[1     2     5     4    10    11    14    13     3     8     6     7    12    17    15    16    19    20    23   22     9    21    26    24    25    18    27];
            %     xi=param_coords(1);
            %     xis = [(xi-1)*xi/2;  (xi+1)*xi/2;  -(xi+1)*(xi-1)];
            %     dxis = [x1; (xi+1/2); -2*xi];
            %     eta=param_coords(2);
            %     etas = [(eta-1)*eta/2;  (eta+1)*eta/2;  -(eta+1)*(eta-1)];
            %     detas=[(eta-1/2); (eta+1/2); -2*eta];
            %     zeta=param_coords(3);
            %     zetas = [(zeta-1)*zeta/2;  (zeta+1)*zeta/2;  -(zeta+1)*(zeta-1)];
            %     dzetas=[(zeta-1/2); (zeta+1/2); -2*zeta];
            %
            %     Nxi =zeros(3, 3, 3);
            %     xisetas=(dxis*etas');
            %     Nxi(:,:,1)=xisetas*zetas(1);
            %     Nxi(:,:,2)=xisetas*zetas(2);
            %     Nxi(:,:,3)=xisetas*zetas(3);
            %
            %     Net =zeros(3, 3, 3);
            %     xisetas=(xis*detas');
            %     Net(:,:,1)=xisetas*zetas(1);
            %     Net(:,:,2)=xisetas*zetas(2);
            %     Net(:,:,3)=xisetas*zetas(3);
            %
            %     Nze =zeros(3, 3, 3);
            %     xisetas=(xis*etas');
            %     Nze(:,:,1)=xisetas*dzetas(1);
            %     Nze(:,:,2)=xisetas*dzetas(2);
            %     Nze(:,:,3)=xisetas*dzetas(3);
            %
            % %     val = [Nxi(idx)',Net(idx)',Nze(idx)'];
            %     val = [Nxi(idx);Net(idx);Nze(idx)]';
            %     return;
            % end
            xi=param_coords(1);
            eta=param_coords(2);
            zet=param_coords(3);
            x1 =(xi-1/2);
            x2 =(xi+1/2);
            x3 =(xi-1);
            x4 =(xi+1);
            z1 =(zet-1);
            z2 =(zet-1/2);
            z3 =(zet+1);
            z4 =(zet+1/2);
            y1 =(eta-1);
            y2 =(eta-1/2);
            y3 =(eta+1);
            y4 = (eta+1/2);
            val = [
                [       1/4*z1*zet*x1*y1*eta,        1/4*z1*zet*x3*xi*y2,        1/4*z2*x3*xi*y1*eta]
                [       1/4*z1*zet*x2*y1*eta,        1/4*z1*zet*x4*xi*y2,        1/4*z2*x4*xi*y1*eta]
                [       1/4*z1*zet*x2*y3*eta,        1/4*z1*zet*x4*xi*y4,        1/4*z2*x4*xi*y3*eta]
                [       1/4*z1*zet*x1*y3*eta,        1/4*z1*zet*x3*xi*y4,        1/4*z2*x3*xi*y3*eta]
                [       1/4*z3*zet*x1*y1*eta,        1/4*z3*zet*x3*xi*y2,        1/4*z4*x3*xi*y1*eta]
                [       1/4*z3*zet*x2*y1*eta,        1/4*z3*zet*x4*xi*y2,        1/4*z4*x4*xi*y1*eta]
                [       1/4*z3*zet*x2*y3*eta,        1/4*z3*zet*x4*xi*y4,        1/4*z4*x4*xi*y3*eta]
                [       1/4*z3*zet*x1*y3*eta,        1/4*z3*zet*x3*xi*y4,        1/4*z4*x3*xi*y3*eta]
                [            -1/2*z1*zet*xi*y1*eta,   1/2*z1*zet*(-x4)*x3*y2,   1/2*z2*(-x4)*x3*y1*eta]
                [  1/2*z1*zet*x2*(-y3)*y1,             -1/2*z1*zet*x4*xi*eta,   1/2*z2*x4*xi*(-y3)*y1]
                [            -1/2*z1*zet*xi*y3*eta,   1/2*z1*zet*(-x4)*x3*y4,   1/2*z2*(-x4)*x3*y3*eta]
                [  1/2*z1*zet*x1*(-y3)*y1,             -1/2*z1*zet*x3*xi*eta,   1/2*z2*x3*xi*(-y3)*y1]
                [            -1/2*z3*zet*xi*y1*eta,   1/2*z3*zet*(-x4)*x3*y2,   1/2*z4*(-x4)*x3*y1*eta]
                [  1/2*z3*zet*x2*(-y3)*y1,             -1/2*z3*zet*x4*xi*eta,   1/2*z4*x4*xi*(-y3)*y1]
                [            -1/2*z3*zet*xi*y3*eta,   1/2*z3*zet*(-x4)*x3*y4,   1/2*z4*(-x4)*x3*y3*eta]
                [  1/2*z3*zet*x1*(-y3)*y1,             -1/2*z3*zet*x3*xi*eta,   1/2*z4*x3*xi*(-y3)*y1]
                [  1/2*(-z3)*z1*x1*y1*eta,   1/2*(-z3)*z1*x3*xi*y2,             -1/2*zet*x3*xi*y1*eta]
                [  1/2*(-z3)*z1*x2*y1*eta,   1/2*(-z3)*z1*x4*xi*y2,             -1/2*zet*x4*xi*y1*eta]
                [  1/2*(-z3)*z1*x2*y3*eta,   1/2*(-z3)*z1*x4*xi*y4,             -1/2*zet*x4*xi*y3*eta]
                [  1/2*(-z3)*z1*x1*y3*eta,   1/2*(-z3)*z1*x3*xi*y4,             -1/2*zet*x3*xi*y3*eta]
                [           -z1*zet*xi*(-y3)*y1,            -z1*zet*(-x4)*x3*eta,  z2*(-x4)*x3*(-y3)*y1]
                [           -(-z3)*z1*xi*y1*eta,  (-z3)*z1*(-x4)*x3*y2,            -zet*(-x4)*x3*y1*eta]
                [ (-z3)*z1*x2*(-y3)*y1,            -(-z3)*z1*x4*xi*eta,            -zet*x4*xi*(-y3)*y1]
                [           -(-z3)*z1*xi*y3*eta,  (-z3)*z1*(-x4)*x3*y4,            -zet*(-x4)*x3*y3*eta]
                [ (-z3)*z1*x1*(-y3)*y1,            -(-z3)*z1*x3*xi*eta,            -zet*x3*xi*(-y3)*y1]
                [           -z3*zet*xi*(-y3)*y1,            -z3*zet*(-x4)*x3*eta,  z4*(-x4)*x3*(-y3)*y1]
                [    -2*(-z3)*z1*xi*(-y3)*y1,     -2*(-z3)*z1*(-x4)*x3*eta,     -2*zet*(-x4)*x3*(-y3)*y1]
                ];
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
            conn = self.conn; % connectivity
            x = gather_values(context.x, conn); % coordinates of nodes
            u = gather_values(context.u, conn); % coordinates of nodes
            % these are all defined in clockwise order, which is what Matlab expects
            ef=[1,9,2,18,6,13,5;
                2,10,3,19,7,14,6;
                3,11,4,20,8,15,7;
                4,12,1,17,5,16,8];
            f=[  1     9    21
                9     2    10
                21    10     3
                12    21    11
                5    16    26
                16     8    15
                26    15     7
                13    26    14
                1    17    22
                17     5    13
                22    13     6
                9    22    18
                2    18    23
                18     6    14
                23    14     7
                10    23    19
                3    19    24
                19     7    15
                24    15     8
                11    24    20
                4    20    25
                20     8    16
                25    16     5
                12    25    17
                21    12     1
                10    21     9
                3    11    21
                11     4    12
                26    13     5
                15    26    16
                7    14    26
                14     6    13
                22     9     1
                13    22    17
                6    18    22
                18     2     9
                23    10     2
                14    23    18
                7    19    23
                19     3    10
                24    11     3
                15    24    19
                8    20    24
                20     4    11
                25    12     4
                16    25    20
                5    17    25
                17     1    12];
            edgecolor = 'Black';
            if (isfield(context,'edgecolor'))
                edgecolor =context.edgecolor;
            end
            if (isfield(context,'colorfield'))
            else
                facecolor = 'red';
                if (isfield(context,'facecolor'))
                    facecolor = context.facecolor;
                end
                edgecolor = 'Black';
                if (isfield(context,'edgecolor'))
                    edgecolor = context.edgecolor;
                end
                options=struct ('facecolor',facecolor,'edgecolor',edgecolor);
            end
            conns = self.conn; % connectivity
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
                    options=struct ('colors',colors);
                end
                if (isfield(context,'edgecolor'))
                    options.edgecolor =context.edgecolor;
                end
                edgecolor='black';
                options.edgecolor =edgecolor;
                if ~strcmp(edgecolor,'none')
                    for i=1:size(ef,1)
                        draw_polyline (gv, xu, conns(:,ef(i,:)), options);
                    end
                end
                options.edgecolor ='none';
                for i=1:size(f,1)
                    draw_polygon (gv, xu, conns(:,f(i,:)), options);
                end
            else
                for ac=1:size(conns,1)
                    conn =conns(ac,:);
                    x = gather_values(context.x, conn); % coordinates of nodes
                    u = gather_values(context.u, conn); % coordinates of nodes
                    if (isfield(context,'colorfield'))
                        colors =gather(context.colorfield, conn, 'values', 'noreshape'); % coordinates of nodes
                        options=struct ('colors',colors);
                    end
                    xu=x+u;
                    if isfield(context,'shrink') & context.shrink~=1
                        pc= [-1,-1,-1;
                            1,-1,-1;
                            1,1,-1;
                            -1,1,-1;-1,-1,1;
                            1,-1,1;
                            1,1,1;
                            -1,1,1;0,-1,-1;
                            1,0,-1;
                            0,1,-1;
                            -1,0,-1;0,-1,1;
                            1,0,1;
                            0,1,1;
                            -1,0,1;-1,-1,0;
                            1,-1,0;
                            1,1,0;
                            -1,1,0;0,0,-1;
                            0,-1,0;
                            1,0,0;
                            0,1,0;
                            -1,0,0;
                            0,0,1;0,0,0];
                        n=27;
                        c=ones(n,1)*(sum(pc)/n);
                        pc=c+context.shrink*(pc-c);
                        xus=xu;
                        for j=1:n
                            N=bfun(self,pc(j,:));
                            xus(j,:)=N'*xu;
                        end
                        xu=xus;
                    end
                    options.edgecolor =edgecolor;
                    if ~strcmp(edgecolor,'none')
                        for i=1:size(ef,1)
                            draw_polyline (gv, xu, ef(i,:), options);
                        end
                    end
                    options.edgecolor ='none';
                    for i=1:size(f,1)
                        draw_polygon (gv, xu, f(i,:), options);
                    end
                end;
            end
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
            [f,v]= isosurf_faces(self, conns, context.x, context.u, context.scalarfield, context.isovalue, 5);
            
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
            val = [conn(:,[1, 4, 3, 2, 12, 11, 10, 9, 21]);
                conn(:,[1, 2, 6, 5, 9, 18, 13, 17, 22]);
                conn(:,[2, 3, 7, 6, 10, 19, 14, 18, 23]);
                conn(:,[3, 4, 8, 7,  11, 20, 15, 19, 24]);
                conn(:,[4, 1, 5, 8, 12, 17, 16, 20, 25]);
                conn(:,[6, 7, 8, 5, 14, 15, 16, 13, 26])];
        end
        
        function val =boundary_fe(self)
            % Get the constructor of the class of the  boundary finite element.
            val = @fe_set_Q9;
        end
        
    end
    
end


