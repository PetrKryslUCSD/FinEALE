% Compute the capacity matrix.
%
% function C = capacity(self, assembler, geom, temp)
%
% Return C as the assembled matrix.

% Arguments
%       self  = heat diffusion model
%       assembler = descendent of the sysmat_assembler class
%       geom=geometry field
%       temp=temperature field
%
function C = capacity(self, assembler, geom, temp)
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
    Cedim =temp.dim*fes.nfens;
    start_assembly(assembler, Cedim, Cedim, size(conns,1), temp.nfreedofs, temp.nfreedofs);
    % Now loop over all fes in the block
    for i=1:size(conns,1)
        conn =conns(i,:);
        x=xs(conn,:);
        Ce =zeros(Cedim);;
        dofnums =reshape (temp,gather_dofnums(temp,conn));
        for j=1:npts
            J = Jacobian_matrix(fes,Nders{j},x);
            Jac = Jacobian_volume(fes,conn, Ns{j}, J, x);
            Ce = Ce + Ns{j}*(Ns{j}' * (cv*Jac * w(j)));%'
        end
        assembler.assemble_symmetric(Ce, dofnums);
    end
    C = make_matrix (assembler);
end
