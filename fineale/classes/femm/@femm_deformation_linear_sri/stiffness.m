function K = stiffness (self, assembler, geom, u)
% Compute the stiffness matrix by computing and assembling the
% matrices of the individual FEs.
%
% function K = stiffness (self, assembler, geom, u)
%
% Return a stiffness matrix.
%     assembler = descendent of the sysmat_assembler class
%     geom=geometry field
%     u=displacement field
    fes = self.fes;% grab the finite elements to work on
    % Integration rule: compute the data needed for  numerical quadrature
    [npts_v Ns_v Nders_v w_v] = integration_data_volumetric (self);
    [npts_s Ns_s Nders_s w_s] = integration_data_deviatoric (self);
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
        D_constrained = tangent_moduli(mat,struct('kind',self.split));
        D_unconstrained = tangent_moduli(mat,struct('kind',[self.split '_shear']));
    end
    % Retrieve data for efficiency
    conns = fes.conn; % connectivity
    labels = fes.label; % connectivity
    xs =geom.values;
    % Prepare assembler
    Kedim =u.dim*fes.nfens;
    start_assembly(assembler, Kedim, Kedim, size(conns,1), u.nfreedofs, u.nfreedofs);
    % Now loop over all fes in the block
    for i=1:size(conns,1)
        conn =conns(i,:);
        x=xs(conn,:);
        dofnums =reshape(u,gather_dofnums(u,conn));
        Ke =zeros(Kedim);
        for j=1:npts_v % loop over all volumetric-strain quadrature points
            c =Ns_v{j}'*x;% physical location of the quadrature point
            J = x' * Nders_v{j};% We compute the Jacobian matrix
            if (~Rm_constant)% do I need to evaluate the local material orientation?
                if (~isempty(labels )),  Rm =Rmh(c,J,labels(i));
                else,                    Rm =Rmh(c,J,[]);                end
            end
            Ndersp = Nders_v{j}/J;% derivatives wrt global coor
            Jac = Jacobian_volume(fes,conn, Ns_v{j}, J, x);
            %if (Jac<=0),error('Non-positive Jacobian');end
            B = self.hBlmat(self,Ns_v{j},Ndersp*Rm,c,Rm);% strain-displacement
            if (~D_constant)
                D_constrained = tangent_moduli(mat,struct('xyz',c,'kind',self.split));
            end
            Ke = Ke + (B'*(D_constrained*(Jac*w_v(j)))*B);
        end
        for j=1:npts_s % loop over all deviatoric-strain quadrature points
            c =Ns_s{j}'*x;% physical location of the quadrature point
            J = x' * Nders_s{j};% We compute the Jacobian matrix
            if (~Rm_constant)% do I need to evaluate the local material orientation?
                if (~isempty(labels )),  Rm =Rmh(c,J,labels(i));
                else,                    Rm =Rmh(c,J,[]);                end
            end
            Ndersp = Nders_s{j}/J;% derivatives wrt global coor
            Jac = Jacobian_volume(fes,conn, Ns_s{j}, J, x);
            %if (Jac<=0),error('Non-positive Jacobian');end
            B = self.hBlmat(self,Ns_s{j},Ndersp*Rm,c,Rm);% strain-displacement
            if (~D_constant)
                D_unconstrained = tangent_moduli(mat,struct('xyz',c,'kind',[self.split '_shear']));
            end
            Ke = Ke + (B'*(D_unconstrained*(Jac*w_s(j)))*B);
        end
        assemble_symmetric(assembler, Ke, dofnums);
    end
    K = make_matrix (assembler);
end
