function [F,self] = restoring_force(self, assembler, geom, un1, un, dt)
% Compute the restoring force vector.
%
%  Call as
%     geom     - geometry field
%     un1      - displacement field at the end of time step t_n+1
%     un       - displacement field at the end of time step t_n
%     dt       - time step from  t_n to t_n+1; needed only by some
%                materials
%
% Note: This method *UPDATES* the state of the FEMM object.  In
%       particular, the material state gets updated.  If this gets
%       called for converged u1, the FEMM must be assigned to itself
%       on return; otherwise the second return value must be ignored!
%
if (~exist('dt','var'))
    % If the time step is needed, this will cause trouble.  Which indicates
    % it should have been supplied. 
    dt=[];
end

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
Uns = un.values; % displacement in step n
Un1s = un1.values; % displacement in step n+1
context.Fn1 = []; context.Fn = []; context.dt=dt;
context.output='Cauchy';
context.get_kinem_var=@get_kinem_var;
% Prepare assembler
Kedim =un1.dim*self.fes.nfens;
start_assembly(assembler, un1.nfreedofs);
for i=1:size(conns,1)
    conn =conns(i,:); % connectivity
    dofnums =reshape(un1,gather_dofnums(un1,conn));
    X=Xs(conn,:);
    Un=Uns(conn,:);
    Un1=Un1s(conn,:);
    xn = X + Un; % previous known coordinates
    xn1 = X + Un1; % current coordinates
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
        Fn1 = xn1'*gradNX;% Current deformation gradient
        context.Fn1=Rm'*Fn1*Rm;%  deformation gradient  in material coordinates
        Fn = xn'*gradNX;% Current deformation gradient
        context.Fn=Rm'*Fn*Rm;%  deformation gradient  in material coordinates
        [cauchy,self.matstates{i,j}] = state(self.material, self.matstates{i,j}, context);
        gcauchy =self.material.stress_vector_rotation(Rm')*cauchy; % material to global
        B = self.hBlmat(self,Ns{j},gradNX/Fn1,c,[]);% strain-displacement
        Fe = Fe - B'* (gcauchy * (Jac * w(j) * det(Fn1))); % note the sign
    end
    assemble(assembler, Fe, dofnums);
end
F = make_vector (assembler);
return
% 
function T=get_kinem_var(type)
    switch (type)
        otherwise
            T=[];
    end
end

end