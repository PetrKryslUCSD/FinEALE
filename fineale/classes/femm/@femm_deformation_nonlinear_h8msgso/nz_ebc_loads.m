function F = nz_ebc_loads(self, assembler, geom, u1)
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

if isempty(self.phis)%  Do we need to calculate the form factors?
    self=self.update(geom, u1, u1);
end

% Integration rule
[npts, Ns, gradNparams, w] = integration_data (self);;
gradN =cell(8,1); Jac=cell(8,1);
% Material orientation
Rm_constant = is_material_orientation_constant(self);% if not constant, need to compute  at each point
if (~Rm_constant)
    Rmh = self.Rm;% handle to a function  to evaluate Rm
else
    Rm = self.Rm;
    if (isempty(Rm)),Rm=eye(3);end % for identity transformation
end
% Retrieve data for efficiency
conns = self.fes.conn; % connectivity
labels = self.fes.label; % connectivity
Xs = geom.values; % reference coordinates
% Us = u.values; % displacement in step n
U1s = u1.values; % displacement in step n+1
context.F = [];
context.stiff_type='Eulerian';
% Prepare assembler
Kedim =u1.dim*self.fes.nfens;
start_assembly(assembler, u1.nfreedofs);
for i=1:size(conns,1)
    conn =conns(i,:); % connectivity
    pu = reshape(u1, gather_fixed_values(u1, conn));
    if norm (pu) ~= 0
        X=Xs(conn,:);
        U1=U1s(conn,:);
        x1 = X + U1; % current coordinates
        dofnums =reshape(u1,gather_dofnums(u1,conn));
        Ke =zeros(Kedim);;
        % First we calculate  the mean basis function gradient matrix and the volume of the element
        gradN_mean =zeros(self.fes.nfens,u1.dim); V=0;
        for j=1:npts % loop over all quadrature points
            J = X' * gradNparams{j};% Jacobian matrix wrt reference coordinates
            gradN{j} = gradNparams{j}/J;% derivatives wrt global reference coor
            Jac{j} = det(J);        if (Jac{j}<=0),error('Non-positive Jacobian');end
            dV=(Jac{j}*w(j));
            gradN_mean =gradN_mean+gradN{j}*dV;
            V =V+dV;
        end
        gradN_mean=gradN_mean/V; % Finish the calculation of the mean basis function gradient matrix
        %     Material orientation?
        if (~Rm_constant)% do I need to evaluate the local material orientation?
            c =mean(X);% physical location of the quadrature point
            if (~isempty(labels )),  Rm =Rmh(c,[],labels(i));%  No Jacobian matrix?
            else,                    Rm =Rmh(c,[],[]);                end
        end
        % Now we calculate the mean deformation gradient
        F1bar =x1' * gradN_mean;
        Bbar = self.hBlmat(self,[],gradN_mean/F1bar,[],[]);% strain-displacement d/dX
        context.ms=self.matstates{i};
        context.F=Rm'*F1bar*Rm;%  Deformation gradient in material coordinates
        D = tangent_moduli (self.material, context);
        D = rotate_stiffness(self.material,D,Rm');
        Ke = (Bbar'*(D*(V))*Bbar);
        context.ms=[];
        Dstab = tangent_moduli (self.stabilization_material, context);
        Dstab = rotate_stiffness(self.stabilization_material,Dstab,Rm');
        Ke = Ke - (Bbar'*(Dstab*(self.phis(i)*V))*Bbar);
        for j=1:npts       % Loop over all integration points
            F1 = x1'*gradN{j};% Current deformation gradient
            B = self.hBlmat(self,[],gradN{j}/F1,[],[]);% strain-displacement
            context.F=Rm'*F1*Rm;
            Dstab = tangent_moduli (self.stabilization_material, context);
            Dstab = rotate_stiffness(self.stabilization_material,Dstab,Rm');
            Ke = Ke + (B'*(Dstab*(self.phis(i)*Jac{j}*w(j)))*B);
        end
        Fe =  -Ke*pu;
        assemble(assembler, Fe, dofnums);
    end
end
F = make_vector (assembler);
end
