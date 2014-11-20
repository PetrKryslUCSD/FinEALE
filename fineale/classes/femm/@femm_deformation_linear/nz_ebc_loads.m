function F = nz_ebc_loads(self, assembler, geom, u)
 % Compute load vector for nonzero essential boundary conditions.
%
% function F = nz_ebc_loads(self, assembler, geom, u)
% 
%
% Return the assembled system vector F.
%    Arguments
%     assembler =  descendent of sysvec_assembler
%     geom=geometry field
%     u=displacement field
%
   fes = self.fes;% grab the finite elements to work on
    % Integration rule: compute the data needed for  numerical quadrature
    [npts Ns Nders w] = integration_data (self);
    % Material orientation
    Rm_constant = is_material_orientation_constant(self);% if not constant, need to compute  at each point
    if (~Rm_constant)
        Rmh = self.Rm;% handle to a function  to evaluate Rm
    else
        Rm = self.Rm;
    end     
    % Material
    mat = self.material;
    D_constant = are_tangent_moduli_constant (mat);
    if (D_constant)
        D = tangent_moduli(mat,[]);
    end
    % Retrieve data for efficiency
    conns = fes.conn; % connectivity
    xs =geom.values;
    labels = fes.label; % connectivity
    % Prepare assembler
    Kedim =u.dim*fes.nfens;
    start_assembly(assembler, u.nfreedofs);
    % Now loop over all fes in the block
    for i=1:size(conns,1)
        conn =conns(i,:);
        pu = reshape(u, gather_fixed_values(u, conn));
        if norm (pu) ~= 0
            x=xs(conn,:);
            dofnums =reshape(u,gather_dofnums(u,conn));
            Ke =zeros(Kedim);;
            for j=1:npts
                c =Ns{j}'*x;% physical location of the quadrature point
                J = x' * Nders{j};% We compute the Jacobian matrix
                if (~Rm_constant)% do I need to evaluate the local material orientation?
                    if (~isempty(labels )),  Rm =Rmh(c,J,labels(i));
                    else,                    Rm =Rmh(c,J,[]);                end
                end
                Ndersp = Nders{j}/J;% derivatives wrt global coor
                Jac = Jacobian_volume(fes,conn, Ns{j}, J, x);
                if (Jac<=0),error('Non-positive Jacobian');end
                B=self.hBlmat(self,Ns{j},Ndersp*Rm,c,Rm);%  strains in material coordinates, displacements in global coordinates
                if (~D_constant)  will% Need to compute the material stiffness at c
                    D = tangent_moduli(mat,struct('xyz',c));%  Moduli in material orientation
                end
                Ke = Ke + (B'*(D*(Jac*w(j)))*B);
            end
            Fe =  -Ke*pu;
            assemble(assembler, Fe, dofnums);
        end    
    end
    F = make_vector (assembler);
end
