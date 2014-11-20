function F = nz_ebc_loads_surface_transfer(self, assembler, geom, temp)
% Compute the load vector produced by nonzero EBC on surface transfer elements.
%
% function F = nz_ebc_loads_surface_transfer(self, assembler, geom, temp)
%
% Compute the load vector produced by nonzero fixed temperature
% at nodes  shared by elements with surface heat transfer condition.
%
% Arguments
%     assembler = descendent of the sysvec_assembler class
%     geom=geometry field
%     temp=temperature field
%
% Return the assembled load vector.
   fes = self.fes;
    % Integration rule
    % Precompute basis f. values + basis f. gradients wrt parametric coor
    [npts Ns Nders w] = integration_data (self);
    % surface transfer coefficient
    h = self.surface_transfer;
    % Prepare some data: 
    conns = fes.conn;% connectivity
    xs =geom.values;% retrieve the geometry information
    % Prepare assembler
    Hedim =temp.dim*fes.nfens;
    start_assembly(assembler, temp.nfreedofs);
    for i=1:size(conns,1)
        conn =conns(i,:);
        pT = reshape(temp, gather_fixed_values(temp, conn));
        if norm (pT) ~= 0
            x=xs(conn,:);
            He =zeros(Hedim);;
            dofnums =reshape (temp,gather_dofnums(temp,conn));
            for j=1:npts
                J = Jacobian_matrix(fes,Nders{j},x);
                Jac = Jacobian_surface(fes,conn, Ns{j}, J, x);
                He = He + ((h*Jac*w(j))*Ns{j})*Ns{j}';%'
            end% Loop over quadrature points
            Fe =  -He*pT;
            assemble(assembler, Fe, dofnums);
        end    
    end% Loop over elements
    F= make_vector (assembler);
end