function H = surface_transfer(self, assembler, geom, temp)
% Compute the surface heat transfer matrix.
%
% function H = surface_transfer(self, assembler, geom, temp)
%
% Arguments
%           self  = heat diffusion model  
%           assembler = descendent of the sysmat_assembler class
%           geom=geometry field
%           temp=temperature field
%
% Returns H as a matrix.
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
    start_assembly(assembler, Hedim, Hedim, size(conns,1), temp.nfreedofs, temp.nfreedofs);
    for i=1:size(conns,1)
        conn =conns(i,:);
        x=xs(conn,:);
        He =zeros(Hedim);;
        dofnums =reshape (temp,gather_dofnums(temp,conn));
        for j=1:npts_per_fe
            J = Jacobian_matrix(fes,Nders{j},x);
            Jac = Jacobian_surface(fes,conn, Ns{j}, J, x);
            He = He + ((h*Jac*w(j))*Ns{j})*Ns{j}';%'
        end
        assembler.assemble_symmetric(He, dofnums);
    end
    H = make_matrix (assembler);
end
