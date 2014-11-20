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
[npts_v Ns_v Nders_v w_v] = integration_data_constrained (self);
[npts_s Ns_s Nders_s w_s] = integration_data_unconstrained (self);
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
    D=tangent_moduli(mat,[]);
    [D_constrained,D_unconstrained]=constrained_mat_data(self,D);
    D_constrained = D_constrained*self.fraction_constrained;
    D_unconstrained = D_unconstrained*self.fraction_unconstrained;
end
% Retrieve data for efficiency
conns = fes.conn; % connectivity
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
        Ke_s =zeros(Kedim);;
        Vfull=0;
        for j=1:npts_s % loop over all deviatoric-strain quadrature points
            c =Ns_s{j}'*x;% physical location of the quadrature point
            J = x' * Nders_s{j};% We compute the Jacobian matrix
            if (~Rm_constant)% do I need to evaluate the local material orientation?
                if (~isempty(labels )),  Rm =Rmh(c,J,labels(i));
                else,                    Rm =Rmh(c,J,[]);                end
            end
            Ndersp = Nders_s{j}/J;% derivatives wrt global coor
            Jac = Jacobian_volume(fes,conn, Ns_s{j}, J, x);
            if (Jac<=0),error('Non-positive Jacobian');end
            if (~D_constant)
                D=tangent_moduli(mat,struct('xyz',c));
                [D_constrained,D_unconstrained]=constrained_mat_data(self,D);
                D_unconstrained = D_unconstrained*self.fraction_unconstrained;
            end
            B=self.hBlmat(self,Ns_s{j},Ndersp*Rm,c,Rm);%  strains in material coordinates, displacements in global coordinates
            Ke_s = Ke_s + (B'*(D_unconstrained*(Jac*w_s(j)))*B);
            Vfull=Vfull+(Jac*w_s(j));
        end
        Ke_v =zeros(Kedim);;
        Vred=0;
        for j=1:npts_v % loop over all volumetric-strain quadrature points
            c =Ns_v{j}'*x;% physical location of the quadrature point
            J = x' * Nders_v{j};% We compute the Jacobian matrix
            if (~Rm_constant)% do I need to evaluate the local material orientation?
                if (~isempty(labels )),  Rm =Rmh(c,J,labels(i));
                else,                    Rm =Rmh(c,J,[]);                end
            end
            Ndersp = Nders_v{j}/J;% derivatives wrt global coor
            Jac = Jacobian_volume(fes,conn, Ns_v{j}, J, x);
            if (Jac<=0),error('Non-positive Jacobian');end
            if (~D_constant)
                D=tangent_moduli(mat,struct('xyz',c));
                [D_constrained,D_unconstrained]=constrained_mat_data(self,D);
                D_constrained = D_constrained*self.fraction_constrained;
           end
            B=self.hBlmat(self,Ns_v{j},Ndersp*Rm,c,Rm);%  strains in material coordinates, displacements in global coordinates
            Ke_v = Ke_v + (B'*(D_constrained*(Jac*w_v(j)))*B);
            Vred=Vred+(Jac*w_v(j));
        end
        Ke=Ke_s+Ke_v*(Vfull/Vred);% account for the insufficiency of the reduced rule  to calculate volume accurately
        Fe =  -Ke*pu;
        assemble(assembler, Fe, dofnums);
    end
end
F = make_vector (assembler);
end
