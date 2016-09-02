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
F=[]; % This is the global vector of restoring forces
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
Uns = un.values; % displacement in step n
Un1s = un1.values; % displacement in step n+1
context.dT = 0;
context.output='Cauchy';
context.dt=dt;% time step
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
    % First we calculate  the mean basis function gradient matrix and the volume of the element
    gradN_mean =zeros(self.fes.nfens,un1.dim); V=0;
    for j=1:npts % loop over all quadrature points
        J = X' * gradNparams{j};% Jacobian matrix wrt reference coordinates
        gradN{j} = gradNparams{j}/J;% derivatives wrt reference coor
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
    % Now we calculate the mean deformation gradient dx/dX   
    Fn1bar = xn1'*gradN_mean;% Current deformation gradient
    context.Fn1=Rm'*Fn1bar*Rm;%  deformation gradient  in material coordinates
    Fnbar = xn'*gradN_mean;% Current deformation gradient
    context.Fn=Rm'*Fnbar*Rm;%  deformation gradient  in material coordinates
    % Update the stress for the real material
    [cauchy,self.matstates{i}] = state(self.material, self.matstates{i}, context);
    % Update the stress for the stabilization material
    [stabcauchy,~] = state(self.stabilization_material, [], context);
    Bbar = self.hBlmat(self,[],gradN_mean/Fn1bar,[],[]);% strain-displacement d/dx
    gcauchy =self.material.stress_vector_rotation(Rm')*(cauchy-self.phis(i)*stabcauchy); % to global
    Fe =  - Bbar'* (gcauchy * (V* det(Fn1bar))) ; % note the sign
    %     Now we update the stress for the stabilization material for the second time; this time for the full quadrature  rule
    for j=1:npts
        Fn1 =xn1'*gradN{j};
        context.Fn1 = Fn1;% Current deformation gradient wrt material orientation
        [stabcauchy,~] = state(self.stabilization_material, [], context);
        B = self.hBlmat(self,[],gradN{j}/Fn1,[],[]);% strain-displacement
        Fe = Fe - B'* (stabcauchy * (self.phis(i)* Jac{j} * w(j)* det(Fn1))) ; % note the sign
    end
    assemble(assembler, Fe, dofnums);
end
F = make_vector (assembler);
return;
end



function [F,self] = restoring_force_for_anisotropic_stabilization_material(self, assembler, geom, un1, un, dt)
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
F=[]; % This is the global vector of restoring forces
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
Uns = un.values; % displacement in step n
Un1s = un1.values; % displacement in step n+1
context.dT = 0;
context.output='Cauchy';
context.dt=dt;% time step
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
    % First we calculate  the mean basis function gradient matrix and the volume of the element
    gradN_mean =zeros(self.fes.nfens,un1.dim); V=0;
    for j=1:npts % loop over all quadrature points
        J = X' * gradNparams{j};% Jacobian matrix wrt reference coordinates
        gradN{j} = gradNparams{j}/J;% derivatives wrt reference coor
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
    % Now we calculate the mean deformation gradient dx/dX   
    Fn1bar = xn1'*gradN_mean;% Current deformation gradient
    context.Fn1=Rm'*Fn1bar*Rm;%  deformation gradient  in material coordinates
    Fnbar = xn'*gradN_mean;% Current deformation gradient
    context.Fn=Rm'*Fnbar*Rm;%  deformation gradient  in material coordinates
    % Update the stress for the real material
    [cauchy,self.matstates{i}] = state(self.material, self.matstates{i}, context);
    % Update the stress for the stabilization material
    [stabcauchy,~] = state(self.stabilization_material, [], context);
    Bbar = self.hBlmat(self,[],gradN_mean/Fn1bar,[],[]);% strain-displacement d/dx
    gcauchy =self.material.stress_vector_rotation(Rm')*(cauchy-self.phis(i)*stabcauchy); % to global
    Fe =  - Bbar'* (gcauchy * (V* det(Fn1bar))) ; % note the sign
    %     Now we update the stress for the stabilization material for the second time; this time for the full quadrature  rule
    for j=1:npts
        Fn1 =xn1'*gradN{j};
        context.Fn1 = Rm'*Fn1*Rm;% Current deformation gradient wrt material orientation
        [stabcauchy,~] = state(self.stabilization_material, [], context);
        B = self.hBlmat(self,[],gradN{j}/context.Fn1,[],[]);% strain-displacement
        gcauchy =self.material.stress_vector_rotation(Rm')*(self.phis(i)*stabcauchy); % to global
        Fe = Fe - B'* (gcauchy * (Jac{j} * w(j)* det(context.Fn1))) ; % note the sign
    end
    assemble(assembler, Fe, dofnums);
end
F = make_vector (assembler);
return;
end

