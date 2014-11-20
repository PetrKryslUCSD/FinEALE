% Finite element (FE) set class.
%
% The set encloses some part of the computational domain.
%
% Class represents a three-manifold finite element set, which encloses some part
% of the computational domain.
% This class is the base class for solid-like finite elements, i.e. every
% usable 3-D finite element has to be derived from it.
%

classdef  fe_set_3_manifold < fe_set
    
    properties (Constant, GetAccess = public)
        dim=3;% manifold dimension (solid)
    end
    
    methods % constructor
        
        function retobj = fe_set_3_manifold(Parameters)
            % Constructor.
            % Parameters:
            % Parameters: same as fe_set.
            if nargin <1
                Parameters = struct([]);
            end
            retobj = retobj@fe_set(Parameters);
        end
        
    end
    
    methods % concrete methods
        
        function Jac = Jacobian_volume(self, conn, N, J, x)
            % Evaluate the volume Jacobian.
            %
            % function Jac = Jacobian_volume(self, conn, N, J, x)
            %
            %   where
            %      self=gcellset,
            %      conn=connectivity of a single cell in which the Jacobian is to be evaluated
            %      N=values of the basis functions, size(N) = [nbfuns,1]
            %      J = Jacobian matrix
            %      x=array of spatial coordinates of the nodes, size(x) = [nbfuns,dim]
            %
            [sdim, ntan] = size(J);
            if     ntan==3
                Jac = det(J);% Compute the Jacobian
            else
                error('Got an incorrect size of tangents');
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
            %      self=gcellset,
            %      conn=connectivity of a single cell in which the Jacobian is to be evaluated
            %      N=values of the basis functions, size(N) = [nbfuns,1]
            %      J = Jacobian matrix
            %      x=array of spatial coordinates of the nodes, size(x) = [nbfuns,dim]
            %
            switch (m)
                case 3
                    Jac = Jacobian_volume(self, conn, N, J, x);
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
            %      self=gcellset,
            %      conn=connectivity of a single cell in which the Jacobian is to be evaluated
            %      N=values of the basis functions, size(N) = [nbfuns,1]
            %      J = Jacobian matrix
            %      x=array of spatial coordinates of the nodes, size(x) = [nbfuns,dim]
            %
            Jac = Jacobian_volume(self, conn, N, J, x);
        end
        
        function draw (self, gv, context)
            % Produce a graphic representation of the three-manifold FE set.
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
            % these are all defined in clockwise order, which is what Matlab expects
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
                context1=context;
                context1.edgecolor ='none';
                for i=1:size(f,1)
                    draw_polygon (gv, xu, conns(:,f(i,:)), context1);
                end
                if (~isfield(context,'edgecolor')) || (~strcmp(context.edgecolor,'none'))
                    for i=1:size(ef,1)
                        draw_polyline (gv, xu, conns(:,ef(i,:)), context);
                    end
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
                    end
                    context1=context;
                    context1.edgecolor ='none';
                    for i=1:size(f,1)
                        draw_polygon (gv, xu, f(i,:), context1);
                    end
                    if (~isfield(context,'edgecolor')) || (~strcmp(context.edgecolor,'none'))
                        for i=1:size(ef,1)
                            draw_polyline (gv, xu, ef(i,:), context);
                        end
                    end
                end
            end
        end
        
        
        % Draw a graphic representation of an isosurface
        % in a three-manifold geometric cell.
        %
        % function draw_isosurface (self, gv, context)
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
        %    scalarfield =field with scalar vertex data,
        %    isovalue = value of the isosurface
        %    color = what color should the isosurface be?  Default: [1,1,1]/2;
        %    nref = how refined should the parametric space of the geometric cell
        %         be? Default: 2.
        %
        function   draw_isosurface(self, gv, context)
            if (isfield(context,'isovalue'))
                isovalue =context.isovalue;
            else
                isovalue =0;
            end
            if ((min(scalars)>isovalue) ||  (max(scalars)<isovalue))
                return
            end
            if (isfield(context,'color'))
                color =context.color;
            else
                color =[1,1,1]/2;
            end
            nref=2;
            if (isfield(context,'nref'))
                nref =context.nref;
            end
            
            xis =linspace(-1,+1,nref);
            [X,Y,Z] = ndgrid(xis,xis,xis);
            V=0*Z;
            for k=1:length(xis)
                for j=1:length(xis)
                    for i=1:length(xis)
                        N = bfun(self, [xis(i),xis(j),xis(k)]);
                        V(i,j,k)=N'*scalars;
                    end
                end
            end
            
            conn = get(self, 'conn'); % connectivity
            if (isempty(conn))
                return; % nothing to be done
            end
            x = gather(context.x, conn, 'values', 'noreshape'); % coordinates of nodes
            u = gather(context.u, conn, 'values', 'noreshape'); % coordinates of nodes
            xu=x+u;
            
            [f,v] = isosurface(X,Y,Z,V,isovalue);
            for i=1:size(v,1)
                N = bfun(self, v(i,:));
                v(i,:)=N'*xu;
            end
            facealpha = 1;
            if (isfield(context,'facealpha'))
                facealpha = context.facealpha;
            end
            patch('Faces',f,'Vertices',v,'FaceColor',color,'EdgeColor','none','FaceAlpha', facealpha);
        end
        
    end
    
end


