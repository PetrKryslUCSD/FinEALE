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
[npts_c,Ns_c,Nders_c,Ns_pv_c,w_c,pc_c] = integration_data_constrained (self);
[npts_u,Ns_u,Nders_u,Ns_pv_u,w_u,pc_u] = integration_data_unconstrained (self);
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
    D = tangent_moduli(mat,[]);% Moduli in material orientation
    [m1,Id,Psi_mat]=constrained_mat_data(self,D);
    srpsi=diag(sqrt(diag(Psi_mat)));
end
PVdim=size(Ns_pv_c{1},1);
Ddim=size(D,1);;
% Retrieve data for efficiency
conns = fes.conn; % connectivity
labels = fes.label; % connectivity
xs =geom.values;

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
        CT =zeros(PVdim,Kedim);
        E=zeros(PVdim,PVdim);
        for j=1:npts_c % loop over all deviatoric-strain quadrature points
            c =Ns_c{j}'*x;% physical location of the quadrature point
            J = x' * Nders_c{j};% We compute the Jacobian matrix
            if (~Rm_constant)% do I need to evaluate the local material orientation?
                if (~isempty(labels )),  Rm =Rmh(c,J,labels(i));
                else,                    Rm =Rmh(c,J,[]);                end
            end
            Ndersp = Nders_c{j}/J;% derivatives wrt global coor
            Jac = Jacobian_volume(fes,conn, Ns_c{j}, J, x);
            if (Jac<=0),error('Non-positive Jacobian');end
            if (~D_constant) 
                D = tangent_moduli(mat,struct('xyz',c));
                [m1,Id,Psi_mat]=constrained_mat_data(self,D);
                srpsi=diag(sqrt(diag(Psi_mat)));
            end
            B = self.hBlmat(self,Ns_c{j},Ndersp*Rm,c,Rm);% strain-displacement
            CT = CT + (Ns_pv_c{j}*(Jac*w_c(j)))*m1'*B;
            E = E + Ns_pv_c{j}*(Ns_pv_c{j}*(Jac*w_c(j)))';
        end
        W=E\CT;
        Ke =zeros(Kedim);;
        for j=1:npts_u % loop over all volumetric-strain quadrature points
            c =Ns_u{j}'*x;% physical location of the quadrature point
            J = x' * Nders_u{j};% We compute the Jacobian matrix
            if (~Rm_constant)% do I need to evaluate the local material orientation?
                if (~isempty(labels )),  Rm =Rmh(c,J,labels(i));
                else,                    Rm =Rmh(c,J,[]);                end
            end
            Ndersp = Nders_u{j}/J;% derivatives wrt global coor
            Jac = Jacobian_volume(fes,conn, Ns_u{j}, J, x);
            if (Jac<=0),error('Non-positive Jacobian');end
            if (~D_constant) 
                D = tangent_moduli(mat,struct('xyz',c));
                [m1,Id,Psi_mat]=constrained_mat_data(self,D);
                srpsi=diag(sqrt(diag(Psi_mat)));
            end
            B = self.hBlmat(self,Ns_u{j},Ndersp*Rm,c,Rm);% strain-displacement
            Bbar= Id*B + m1*(srpsi/3)*Ns_pv_u{j}'*W;
            Ke = Ke + (Bbar'*(D*(Jac*w_u(j)))*Bbar);
        end
        Fe =  -Ke*pu;
        assemble(assembler, Fe, dofnums);
    end
end
F = make_vector (assembler);
end
