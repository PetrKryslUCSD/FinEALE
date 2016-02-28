%  Draw graphic representation for material states in all gcells.
%
% Input arguments
% self = self
% gv = graphic viewer
% context = struct
% with typically mandatory fields
%    x=reference geometry field
%    u=displacement field
% with optional fields
%    scale = scale the representation of the material state graphic by this
%       factor
%    dT=temperature difference (with respect to the reference temperature);
%       temperature differences assumed to be zero if it is not supplied.
% Any other fields in the context struct may be interpreted by the
% implementation of the material: see the draw() function for the
% particular material.
%
function draw_integration_points(self, gv, context)
    error (' Not implemented')
    gcells = get(self.feblock,'gcells');
    ngcells = count(gcells);
    nfens = get(gcells(1),'nfens');
    dim = get(context.x,'dim');
    % Integration rule
    integration_rule = get(self.feblock, 'integration_rule');
    pc = get(integration_rule, 'param_coords');
    w  = get(integration_rule, 'weights');
    npts_per_gcell = get(integration_rule, 'npts'); % number of integration point
    % Material
    mat = get(self.feblock, 'mater');
    matstates = get(self.feblock, 'matstates');
    % Now loop over all gcells in the block
    context.update_context= [];
    for i=1:ngcells
        conn = get(gcells(i), 'conn'); % connectivity
        x = gather(context.x, conn, 'values', 'noreshape'); % reference coordinates
        U = gather(context.u, conn, 'values', 'noreshape'); % displacement
        if isfield(context,'dT')
            dTs = gather(context.dT, conn, 'values', 'noreshape'); % displacement
        else
            dTs = zeros(nfens,1);
        end
        % Loop over all integration points
        for j=1:npts_per_gcell
            Nder = Ndermat_param (gcells(i), pc(j,:));
            ts=Nder'*x; % rows: tangent vectors to parametric curves
            es=local_basis(self,ts); % calculate the local Cartesian basis
            xl=x*es(:,1:2);% local coordinates
            ul=U*es(:,1:2);% local displacements
            [Nspatialder,detF] = Ndermat_spatial (gcells(i), Nder, xl);
            [Nderxl,J] = Ndermat_spatial (gcells(i), Nder, xl);
            xyz = map_to_xyz(gcells(i),pc(j,:),x+U);
            context.ms=matstates{i,j};
            context.xyz=xyz;
            context.es = es;
            context.update_context.gradU = def_grad(gcells(i), Nderxl, ul);
            N = Nmat (gcells(i), pc(j,:));
            context.update_context.dT = N'*dTs;
            draw(mat, gv, context);
        end
    end
end
