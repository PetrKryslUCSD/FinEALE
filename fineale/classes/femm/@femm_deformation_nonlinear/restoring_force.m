function [F,self] = restoring_force(self, assembler, geom, u1, u)
% Compute the restoring force vector.
%
%  Call as
%     geom     - geometry field
%     u1       - displacement field at the end of time step t_n+1
%     u        - displacement field at the end of time step t_n
%
% Note: This method *UPDATES* the state of the FEMM object.  In
%       particular, the material state gets updated.  If this gets
%       called for converged u1, the FEMM must be assigned to itself
%       on return; otherwise the second return value must be ignored!
%
F=[]; % This is the global vector of restoring forces
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
% Us = u.values; % displacement in step n
U1s = u1.values; % displacement in step n+1
context.F = [];
context.output='Cauchy';
% Prepare assembler
Kedim =u.dim*self.fes.nfens;
start_assembly(assembler, u.nfreedofs);
for i=1:size(conns,1)
    conn =conns(i,:); % connectivity
    dofnums =reshape(u1,gather_dofnums(u1,conn));
    X=Xs(conn,:);
%     U=Us(conn,:);
    U1=U1s(conn,:);
    %     x = X + U; % previous known coordinates
    x1 = X + U1; % current coordinates
    context.dT = 0;
    Fe=zeros(Kedim,1);
    for j=1:npts
        c =Ns{j}'*X;% physical location of the quadrature point
        J = X' * gradNs{j};% We compute the Jacobian matrix
        if (~Rm_constant)% do I need to evaluate the local material orientation?
            if (~isempty(labels )),  Rm =Rmh(c,J,labels(i));
            else,                    Rm =Rmh(c,J,[]);                end
        end
        gradNX = gradNs{j}/J;% derivatives wrt global coor
        Jac = Jacobian_volume(self.fes, conn, Ns{j}, J, X);
        if (Jac<=0),error('Non-positive Jacobian');end
        F1 = X'*gradNX + U1'*gradNX;% Current deformation gradient
        context.F=Rm'*F1*Rm;%  deformation gradient  in material coordinates
        [cauchy,self.matstates{i,j}] = update(self.material, self.matstates{i,j}, context);
        gcauchy =self.material.stress_vector_rotation(Rm')*cauchy; % material to global
        B = self.hBlmat(self,Ns{j},gradNX/F1,c,[]);% strain-displacement
        Fe = Fe - B'* (gcauchy * (Jac * w(j) * det(F1))); % note the sign
    end
    assemble(assembler, Fe, dofnums);
end
F = make_vector (assembler);
end