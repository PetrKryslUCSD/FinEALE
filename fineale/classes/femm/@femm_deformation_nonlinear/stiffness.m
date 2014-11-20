function K = stiffness (self, assembler, geom, u1, u)
% Compute the stiffness matrix by computing and assembling the
% matrices of the individual FEs.
%
% function K = stiffness (self, assembler, geom, u1, u)
%
% Return a stiffness matrix.
%     assembler = descendent of the sysmat_assembler class
%     geom=geometry field
%     u1=displacement field, the current configuration
%     u=displacement field, the last known configuration

% Integration rule
[npts Ns gradNs w] = integration_data (self);;
% Material orientation
Rm_constant = is_material_orientation_constant(self);% if not constant, need to compute  at each point
if (~Rm_constant)
    Rmh = self.Rm;% handle to a function  to evaluate Rm
else
    Rm = self.Rm;
end
% Retrieve data for efficiency
conns = self.fes.conn; % connectivity
labels = self.fes.label; % connectivity
Xs = geom.values; % reference coordinates
Us = u.values; % displacement in step n
U1s = u1.values; % displacement in step n+1
context.F = [];
context.stiff_type='Eulerian';
% Prepare assembler
Kedim =u.dim*self.fes.nfens;
start_assembly(assembler, Kedim, Kedim, size(conns,1), u.nfreedofs, u.nfreedofs);
for i=1:size(conns,1)
    conn =conns(i,:); % connectivity
    X=Xs(conn,:);
    U=Us(conn,:);
    U1=U1s(conn,:);
    x = X + U; % previous known coordinates
    x1 = X + U1; % current coordinates
    dofnums =reshape(u1,gather_dofnums(u1,conn));
    Ke =zeros(Kedim);;
    for j=1:npts       % Loop over all integration points
        c =Ns{j}'*X;% physical location of the quadrature point
        J = X' * gradNs{j};% We compute the Jacobian matrix
        if (~Rm_constant)% do I need to evaluate the local material orientation?
            if (~isempty(labels )),  Rm =Rmh(c,J,labels(i));
            else,                    Rm =Rmh(c,J,[]);                end
        end
        gradNX = gradNs{j}/J;% derivatives wrt global coor
        Jac = Jacobian_volume(self.fes, conn, Ns{j}, J, X);
        if (Jac<=0),error('Non-positive Jacobian');end
        F1 = x1'*gradNX;% Current deformation gradientwill
        B = self.hBlmat(self,Ns{j},gradNX/F1,c,[]);% strain-displacement
        context.ms =self.matstates{i,j};
        context.F=Rm'*F1*Rm;%  deformation gradient  in material coordinates
        D = tangent_moduli (self.material, context);% in material coordinates
        D = rotate_stiffness(self.material,D,Rm');%in global coordinates
        Ke = Ke + (B'*(D*(Jac*w(j)*det(F1)))*B);
    end
    assemble_symmetric(assembler, Ke, dofnums);
end
K = make_matrix (assembler);
end
