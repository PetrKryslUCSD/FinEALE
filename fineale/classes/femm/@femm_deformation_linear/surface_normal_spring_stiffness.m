function H = surface_normal_spring_stiffness(self, assembler, geom, u)
% Compute the stiffness matrix of surface normal spring.
%
% function H = surface_normal_spring_stiffness(self, assembler, geom, u)
%
% Arguments
%           self  = heat diffusion model  
%           assembler = descendent of the sysmat_assembler class
%           geom=geometry field
%           u=displacement field
%
% Returns H as a matrix.
%
% Rationale: consider continuously distributed springs between the surface of 
% the solid body and the 'ground', in the direction normal to the surface. 
% If the spring coefficient becomes large, we have an approximate
% method of enforcing the normal displacement to the surface.
    fes = self.fes;
    % Integration rule: compute the data needed for  numerical quadrature
    [npts Ns Nders w] = integration_data (self);
    % surface normal spring coefficient
    k = self.surface_normal_spring_coefficient;
    % Prepare some data: 
    conns = fes.conn;% connectivity
    xs =geom.values;% retrieve the geometry information
    % Prepare assembler
    Hedim =u.dim*fes.nfens;
    start_assembly(assembler, Hedim, Hedim, size(conns,1), u.nfreedofs, u.nfreedofs);
    for i=1:size(conns,1)
        conn =conns(i,:);
        x=xs(conn,:);
        He =zeros(Hedim);;
        dofnums =reshape (u,gather_dofnums(u,conn));
        for j=1:npts
            Y=transpose(Ns{j})*x; % Point at which pressure force acts
            J = Jacobian_matrix(fes,Nders{j},x);
            Jac = Jacobian_surface(fes,conn, Ns{j}, J, x);
            n= normal(self,Y,J);% find the normal to the surface
            Nn =reshape(n*Ns{j}',Hedim,1);% The normal n is a column vector
            He = He + Nn*(Nn'*k*(Jac*w(j)));
        end
        assembler.assemble_symmetric(He, dofnums);
    end
    H = make_matrix (assembler);
end
