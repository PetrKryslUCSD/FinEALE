function draw_integration_points(self, gv, context)
%  Draw graphic representation for material states in all finite elements.
%
% function draw_integration_points(self, gv, context)
%
% Input arguments
% self = self
% gv = graphic viewer
% context = struct
% with typically mandatory fields
%    x=reference geometry field
%    u=displacement field to be used for the calculation of the stresses
%    dT=temperature difference field
% with optional fields
%    scale = scale the representation of the material state graphic by this
%       factor
%    u_scale = scaling of the displacement field to be used for display
%    dT=temperature difference (with respect to the reference temperature);
%       temperature differences assumed to be zero if it is not supplied.
% Any other fields in the context struct may be interpreted by the
% implementation of the material: see the draw() function for the
% particular material.
%
output='Cauchy';
if (isfield(context,'output'))
    output= context.output;
end
scale = 1.0;
if (isfield(context,'scale'))
    scale= context.scale;
end
u_scale = 0.0;
if (isfield(context,'u_scale'))
    u_scale= context.u_scale;
end
map_colors = false;
if isfield(context,'data_cmap')
    map_colors = true;
end
if ~isfield(context,'cheap_arrow')
    context.cheap_arrow = true;
end
    context.ms=[];
    context.F=[];
        context.dT=[];

    function idat=inspector(idat, out, xyz, u, pc)
        switch (output)
            case {'stress','Cauchy'}
                sigmav= out;
                [V,D]=eig(stress_6v_to_3x3t(self.material,sigmav));
                [princStress,ix] =sort(diag(D),'descend');
                V=V(:,ix);
                component =1;
                if isfield(context,'component')
                    component = context.component;
                end
                if map_colors
                    context.color = map_data(context.data_cmap, sigmav(component));
                end
                xyz=(xyz+u_scale*u)';
                if (context.cheap_arrow)
                    a=scale*V(:,1)*princStress(1); if ~map_colors, context.color= 'c'; end
                    draw_arrow(gv, xyz-a, 2*a, context);
                    a=scale*V(:,2)*princStress(2); if ~map_colors, context.color= 'y'; end
                    draw_arrow(gv, xyz-a, 2*a, context);
                    a=scale*V(:,3)*princStress(3); if ~map_colors, context.color= 'm'; end
                    draw_arrow(gv, xyz-a, 2*a, context);
                else
                    a=scale*V(:,1)*princStress(1); if ~map_colors, context.color= 'c'; end
                    if (princStress(1)>0)
                        draw_arrow(gv, xyz, a, context); draw_arrow(gv, xyz, -a, context);
                    else
                        draw_arrow(gv, xyz-a, a, context); draw_arrow(gv, xyz+a, -a, context);
                    end
                    a=scale*V(:,2)*princStress(2); if ~map_colors, context.color= 'y'; end
                    if (princStress(2)>0)
                        draw_arrow(gv, xyz, a, context); draw_arrow(gv, xyz, -a, context);
                    else
                        draw_arrow(gv, xyz-a, a, context); draw_arrow(gv, xyz+a, -a, context);
                    end
                    a=scale*V(:,3)*princStress(3); if ~map_colors, context.color= 'm'; end
                    if (princStress(3)>0)
                        draw_arrow(gv, xyz, a, context); draw_arrow(gv, xyz, -a, context);
                    else
                        draw_arrow(gv, xyz-a, a, context); draw_arrow(gv, xyz+a, -a, context);
                    end
                end
                %                 context.facecolor= color;
                %                 context.tessel=12;
                %                 draw_ellipsoid(gv, xyz, V, scale*diag(D), context);
            case {'equiv_pl_def'}
                color =[1 1 0];
                if isfield(context,'data_cmap')
                    color = map_data(context.data_cmap, out(1));
                end
                context.facecolor= color;
                draw_ellipsoid(gv, xyz+u_scale*u, eye(3), scale*[1,1,1]*out(1), context);
            otherwise
                % Do nothing
        end
    end
    idat= [];
    idat = inspect_integration_points(self, ...
                    context.x, context.un1, context.un, context.dt, context.dT, 1:count(self.fes), context,...
                    @inspector, idat);
end
%
%
%     fes = get(self.femmlock,'gcells');
%     ngcells = count(gcells);
%     nfens = get(gcells(1),'nfens');
%     dim = get(context.x,'dim');
%     % Integration rule
%     integration_rule = get(self.femmlock, 'integration_rule');
%     pc = get(integration_rule, 'param_coords');
%     w  = get(integration_rule, 'weights');
%     npts_per_gcell = get(integration_rule, 'npts'); % number of integration point
%     for j=1:npts_per_gcell
%         Ns{j} = bfun(gcells,pc(j,:));
%         Nders{j} = bfundpar(gcells,pc(j,:));
%     end
%     % Material orientation /directions
%     Rm = get(self,'Rm');
%     eval_Rm = 0;
%     if strcmp(class(Rm),'function_handle')
%         eval_Rm = true; Rmh =Rm;
%     end
%     % Material
%     mat = get(self.femmlock, 'mater');
%     % Now loop over all gcells in the block
%     conns = get(gcells, 'conn'); % connectivity
%     context.update_context= [];
%     context.ms = [];
%     xs =gather(context.x,(1:get(context.x,'nfens')),'values','noreshape');
%     for i=1:ngcells
%         conn =conns(i,:);
%         X = xs(conn,:);% reference coordinates
%         U_displ = gather(context.u_displ, conn, 'values', 'noreshape'); % displacement
%         Uv = gather(context.u, conn, 'values'); % displacement
%         if isfield(context,'dT')
%             dTs = gather(context.dT, conn, 'values', 'noreshape'); % displacement
%         else
%             dTs = zeros(nfens,1);
%         end
%         % Loop over all integration points
%         for j=1:npts_per_gcell
%             c =Ns{j}'*X;
%             J = X' * Nders{j};%J = Jacobian_matrix(gcells,Nders{j},x);%J = x' * Nder;
%             if (eval_Rm)% do I need to evaluate the local material orientation?
%                 Rm =Rmh(c,J);
%             end
%             Jac = Jacobian_volume(gcells,conn, Ns{j}, J, X);
%             if (isempty(Rm))
%                 Ndersp = Nders{j}/J;% derivatives wrt local coor
%             else
%                 Ndersp = Nders{j}/(Rm'*J);% derivatives wrt local coor
%             end
%             B = self.hBlmat(self,Ns{j},Ndersp,c,Rm);% strain-displacement
%             context.D = tangent_moduli(mat,struct('xyz',c));
%             context.update_context.strain = B*Uv;
%             context.update_context.dT =Ns{j}'*dTs;
%             context.xyz =Ns{j}'*(X+U_displ);
%             context.Rm=Rm;
%             draw(mat, gv, context);
%         end
%     end
%     return;
% end
