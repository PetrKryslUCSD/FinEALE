function F = nz_ebc_loads_acceleration(self, assembler, geom, a)
 % Compute load vector for nonzero acceleration essential boundary conditions.
%
% function F = nz_ebc_loads_acceleration(self, assembler, geom, a)
% 
%
% Return the assembled system vector F.
%    Arguments
%     assembler =  descendent of sysvec_assembler
%     geom=geometry field
%     a=acceleration field
%
   fes = self.fes;% grab the finite elements to work on
    % Integration rule: compute the data needed for  numerical quadrature
    [npts Ns Nders w] = integration_data (self);
    % Material
    mat = self.material;
    rho = mat.property.rho;
    nfens=fes.nfens; dim=a.dim;
    Nexp = zeros(a.dim,a.dim*fes.nfens);
    % Retrieve data for efficiency
    conns = fes.conn; % connectivity
    xs =geom.values;
    % Prepare assembler
    Medim =a.dim*fes.nfens;
    start_assembly(assembler, a.nfreedofs);
    for i=1:size(conns,1)
        conn =conns(i,:);
        pa = reshape(a, gather_fixed_values(a, conn));
        if norm (pa) ~= 0
            x=xs(conn,:);
            dofnums =reshape(a,gather_dofnums(a,conn));
            Me =zeros(Medim);
            for j=1:npts
                J = Jacobian_matrix(fes,Nders{j},x);
                Jac = Jacobian_volume(fes,conn, Ns{j}, J, x);
                for l = 1:nfens
                    Nexp(1:dim,(l-1)*dim+1:(l)*dim)=eye(dim)*Ns{j}(l);
                end;
                Me = Me + (Nexp'*Nexp) * (rho * Jac * w(j));
            end
            Fe =  -Me*pa;
            assemble(assembler, Fe, dofnums);
        end
    end
    F = make_vector (assembler);
end