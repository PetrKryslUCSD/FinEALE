function F = nz_ebc_loads_capacity (self, assembler, geom, temp_rate)
% Compute load vectors for nonzero EBC of fixed temperature rate.
%
% function F = nz_ebc_loads_capacity (self, assembler, geom, temp_rate)
%
% Compute the load vector corresponding to nonzero
% essential boundary conditions (fixed temperature rate).
% 
%    Arguments
%     assembler =  descendent of sysvec_assembler
%     geom=geometry field
%     temp=temperature rate field
%
    fes = self.fes;
    % Precompute basis f. values + basis f. gradients wrt parametric coor
    [npts Ns Nders w] = integration_data (self);
    % Material
    mat = self.material;
    cv = mat.property.specific_heat;
    % Prepare some data: 
    conns = fes.conn; % connectivity
    xs = geom.values;
    % Prepare assembler
    Cedim =temp_rate.dim*fes.nfens;
    start_assembly(assembler, temp_rate.nfreedofs);
    % Now loop over all fes in the block
    for i=1:size(conns,1)
        conn =conns(i,:);
        pT = reshape(temp_rate, gather_fixed_values(temp_rate, conn));
        if norm (pT) ~= 0
            x=xs(conn,:);
            Ce =zeros(Cedim);;
            dofnums =reshape(temp_rate,gather_dofnums(temp_rate,conn));
            for j=1:npts
                J = Jacobian_matrix(fes,Nders{j},x);
                Jac = Jacobian_volume(fes,conn, Ns{j}, J, x);
                Ce = Ce + Ns{j}*(Ns{j}' * (cv*Jac * w(j)));%'
            end
            Fe =  -Ce*pT;
            assemble(assembler, Fe, dofnums);
        end
    end
    F= make_vector (assembler);
end
