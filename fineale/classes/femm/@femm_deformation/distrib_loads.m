function F = distrib_loads(self, assembler, geom, u, fi, m)
% Compute the load vector due to distributed force.
%
% function F = distrib_loads(self, assembler, geom, u, fi, m)
%
% Compute the load vector corresponding to applied distributed
% load. Here it means force per unit volume where volume could be
% length^3, length^2, length^1, or length^0, depending on the manifold
% dimension of the finite element.
%
% Return the assembled system vector F.
%
%    Arguments
%     assembler =  descendent of sysvec_assembler
%     geom=geometry field
%     u=displacement field
%     fi=force intensity object
%     m= manifold dimension
%
    fes = self.fes;
    % Integration rule
    [npts Ns Nders w] = integration_data (self);;
    % Retrieve data for efficiency
    conns = fes.conn; % connectivity
    xs =geom.values;
    % Prepare assembler
    Kedim =u.dim*fes.nfens;
    start_assembly(assembler, u.nfreedofs);
    % Now loop over all fes in the block
    for i=1:size(conns,1)
        conn =conns(i,:);
        x=xs(conn,:);
        dofnums =reshape(u,gather_dofnums(u,conn));
        Fe =zeros(Kedim,1);
        for j=1:npts
            J = x' * Nders{j};% We compute the Jacobian matrix
            Jac = Jacobian_mdim(fes,conn, Ns{j}, J, x, m);
            f=get_magn(fi,Ns{j}'*x,J);%'
            Fe = Fe + reshape(f*Ns{j}',Kedim,1) * (Jac * w(j));%'
        end
        assemble(assembler, Fe, dofnums);
    end
    F = make_vector (assembler);
end
