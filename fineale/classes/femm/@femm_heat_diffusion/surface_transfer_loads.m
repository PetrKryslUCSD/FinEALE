function F = surface_transfer_loads (self, assembler, geom, temp, amb)
% Compute the load vector corresponding to surface heat transfer.
%
% function F = surface_transfer_loads (self, assembler, geom, temp, amb)
%
%    Arguments
%     assembler =  descendent of sysvec_assembler
%     geom=geometry field
%     temp=temperature field
%     amb = ambient temperature field
%
% Return the assembled system vector F.
    fes = self.fes;
    % Integration rule
    integration_rule = self.integration_rule;
    pc = integration_rule.param_coords;
    w  =  integration_rule.weights ;
    npts_per_fe = integration_rule.npts;
    % Precompute basis f. values + basis f. gradients wrt parametric coor
    for j=1:npts_per_fe
        Ns{j} = bfun(fes,pc(j,:));
        Nders{j} = bfundpar(fes,pc(j,:));
    end
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
        pT = reshape(temp, gather_fixed_values(amb, conn));
        if norm (pT) ~= 0
            x=xs(conn,:);
            He =zeros(Hedim);;
            dofnums =reshape (temp,gather_dofnums(temp,conn));
            for j=1:npts_per_fe
                J = Jacobian_matrix(fes,Nders{j},x);
                Jac = Jacobian_surface(fes,conn, Ns{j}, J, x);
                He = He + ((h*Jac*w(j))*Ns{j})*Ns{j}';%'
            end% Loop over quadrature points
            Fe =  He*pT;
            assemble(assembler, Fe, dofnums);
        end    
    end% Loop over elements
    F= make_vector (assembler);
end
   