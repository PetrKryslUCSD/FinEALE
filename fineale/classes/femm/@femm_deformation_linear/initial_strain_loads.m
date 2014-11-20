% Compute the initial-strain load vectors of the individual gcells.
%
% function evs = initial_strain_loads(self, geom, u)
%
% Return an array of the element vectors so they may be assembled.
%    Call as
% ems = initial_strain_loads(femm, geom, u)
%     geom=geometry field
%     u=displacement field
%
function evs = initial_strain_loads(self, geom, u)
    gcells = get(self,'gcells');
    ngcells = count(gcells);
    nfens = get(gcells(1),'nfens');
    dim = get(geom,'dim');
    conns = get(gcells, 'conn'); % connectivity
    % Integration rule
    integration_rule = get(self.femmlock, 'integration_rule');
    pc = get(integration_rule, 'param_coords');
    w  = get(integration_rule, 'weights');
    npts_per_gcell = get(integration_rule, 'npts');
    for j=1:npts_per_gcell
        Ns{j} = bfun(gcells,pc(j,:));
        Nders{j} = bfundpar(gcells,pc(j,:));
    end
    % Material
    mat = get(self, 'mater');
    matstates = get(self, 'matstates');
    %     Preallocate
    Fe = cell(size(conns,1),1);
    eqnums = cell(size(conns,1),1);
    for i=1:size(conns,1)
        eqnums{i} =gather(u,conns(i,:),'eqnums');
        Fe{i} =zeros(get(geom,'dim')*nfens,1);
    end
    xs =gather(geom,(1:get(geom,'nfens')),'values','noreshape');
    % Now loop over all gcells in the block
    for i=1:ngcells
        conn =conns(i,:);
        x=xs(conn,:);
        U = gather(u, conn, 'values'); % displacements
        % Loop over all integration points
        for j=1:npts_per_gcell
            c =Ns{j}'*x;
            J = Jacobian_matrix(gcells,Nders{j},x);
            Jac = Jacobian_volume(gcells,conn, Ns{j}, J, x);
            Ndersp = Nders{j}/J;
            B = self.hBlmat(self,Ndersp,c,[]);%Rm(i,:,:));
            context.D = tangent_moduli(mat,struct('xyz',c));
            context.ms =matstates{i,j};
            iSigma = initial_stress(mat,context);
            Fe{i} = Fe{i} - B'*(iSigma * (Jac * w(j)));
        end
    end
    evs  = elevecset(struct('vec',{Fe}, 'eqnums',{eqnums}));
end
