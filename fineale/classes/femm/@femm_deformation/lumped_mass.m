function M = lumped_mass (self, assembler, geom, u)
% Compute the lumped (diagonal) mass matrix.
%
% function M = lumped_mass (self, assembler, geom, u)
%
% Return a lumped (diagonal) mass matrix.
%     assembler = descendent of the sysmat_assembler class
%     geom=geometry field
%     u=displacement field
%
    fes = self.fes;% grab the finite elements to work on
    % Integration rule: compute the data needed for  numerical quadrature
    [npts Ns Nders w] = integration_data (self);
    % Material
    mat = self.material;
    rho = mat.property.rho;
    nfens=fes.nfens; dim=u.dim;
    Nexp = zeros(u.dim,u.dim*fes.nfens);
    % Retrieve data for efficiency
    conns = fes.conn; % connectivity
    xs =geom.values;
    % Prepare assembler
    Medim =u.dim*fes.nfens;
    start_assembly(assembler, Medim, Medim, size(conns,1), u.nfreedofs, u.nfreedofs);
    for i=1:size(conns,1)
        conn =conns(i,:);
        x=xs(conn,:);
        eqnums =reshape(u,gather_dofnums(u,conn));
        Me =zeros(Medim);
        for j=1:npts
            J = Jacobian_matrix(fes,Nders{j},x);
            Jac = Jacobian_volume(fes,conn, Ns{j}, J, x);
            for l = 1:nfens
                Nexp(1:dim,(l-1)*dim+1:(l)*dim)=eye(dim)*Ns{j}(l);
            end;
            Me = Me + (Nexp'*Nexp) * (rho * Jac * w(j));
        end
        % Hinton,
        em2=sum(sum(Me));
        dem2=sum(diag(Me));
        Me=diag(diag(Me)/dem2*em2);
        assemble_symmetric(assembler, Me, eqnums);
    end
    M = make_matrix (assembler);
end
