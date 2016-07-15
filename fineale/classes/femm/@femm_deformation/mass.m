function M = mass (self, assembler, geom, u)
% Compute the mass matrix.
%
% function M = mass (self, assembler, geom, u)
%
% Return a mass matrix.
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
    % Retrieve data for efficiency
    conns = fes.conn; % connectivity
    xs =geom.values;
    %     Precompute for efficiency
    Nexp = zeros(u.dim,u.dim*fes.nfens);
    NexpNexp={};
    for j=1:npts
        for l = 1:nfens
            Nexp(1:dim,(l-1)*dim+1:(l)*dim)=eye(dim)*Ns{j}(l);
        end;
        NexpNexp{j}=Nexp'*Nexp;
    end
    % Prepare assembler
    Medim =u.dim*fes.nfens;
    start_assembly(assembler, Medim, Medim, size(conns,1), u.nfreedofs, u.nfreedofs);
    for i=1:size(conns,1)
        conn =conns(i,:);
        x=xs(conn,:);
        dofnums =reshape(u,gather_dofnums(u,conn));
        Me =zeros(Medim);
        for j=1:npts
            J = Jacobian_matrix(fes,Nders{j},x);
            Jac = Jacobian_volume(fes,conn, Ns{j}, J, x);
            Me = Me + (NexpNexp{j}) * (rho * Jac * w(j));
        end
        assemble_symmetric(assembler, Me, dofnums);
    end
    M = make_matrix (assembler);
end
