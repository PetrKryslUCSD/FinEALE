function F = distrib_loads(self, assembler, geom, temp, fi, m)
% Compute the distributed-load vector.
%
% function F = distrib_loads(self, assembler, geom, temp, fi, m)
%
% Return the assembled vector due to either internal heat generation 
% (load per unit volume), or due to applied heat flux on the surface.
%
% Arguments
%       assembler =  descendent of sysvec_assembler
%       geom=geometry field
%       temp=temperature field
%       fi=force intensity object
%       m= manifold dimension, 2= surface, 3= volume
    fes = self.fes;
    % Precompute basis f. values + basis f. gradients wrt parametric coor
    [npts Ns Nders w] = integration_data (self);
    % Material
    mat = self.material;
    start_assembly(assembler, temp.nfreedofs);
    % Now loop over all fes in the block
    conns = fes.conn; % connectivity
    xs =geom.values;
    for i=1:size(conns,1)
        conn =conns(i,:);
        x=xs(conn,:);
        dofnums =reshape(temp,gather_dofnums(temp,conn));
        Fe =zeros(fes.nfens,1);
        for j=1:npts
            J = Jacobian_matrix(fes,Nders{j},x);
            Jac = Jacobian_mdim(fes,conn, Ns{j}, J, x, m);
            f=get_magn(fi,transpose(Ns{j})*x,J);
            Fe = Fe + Ns{j} *  (f * Jac * w(j));
        end
        assemble(assembler, Fe, dofnums);
    end
    F= make_vector (assembler);
end
