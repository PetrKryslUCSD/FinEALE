function [self] = update(self,  geom, u1, u)
% Update the state of the FEMM.
%
%  Call as
%     geom     - geometry field
%     u1       - displacement field at the end of time step t_n+1
%     u        - displacement field at the end of time step t_n
%
% Note: This method *UPDATES* the state of the FEMM object.  In
%       particular, the material state gets updated.  This method is
%       supposed to be called for converged u1, the FEMM must be assigned to itself
%       on return.
%
% Note: The state of the FEMM may depend on the geometry field.  It is
%       assumed that  until  the next update() the geometry field geom DOES NOT
%       CHANGE.



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
context.dT = 0;
context.F = [];
context.output='Cauchy';
    % Prepare assembler
    % Kedim =u.dim*self.fes.nfens;
    % start_assembly(assembler, u.nfreedofs);
    
% Compute  +Store the stabilization parameters
self.phis=zeros(size(conns,1),1);;
for i=1:size(conns,1)
    conn =conns(i,:); % connectivity
    dofnums =reshape(u1,gather_dofnums(u1,conn));
    X=Xs(conn,:);
    %     U=Us(conn,:);
    U1=U1s(conn,:);
    %     x = X + U; % previous known coordinates
    x1 = X + U1; % current coordinates
    context.dT = 0;
    %     Fe=zeros(Kedim,1);
    % First we calculate  the mean basis function gradient matrix and the volume of the element
    gradN_mean =zeros(self.fes.nfens,u.dim); V=0;
    self.phis(i)=0;%  stabilization fraction
    for j=1:npts % loop over all quadrature points
        J = X' * gradNparams{j};% Jacobian matrix wrt reference coordinates
        gradN{j} = gradNparams{j}/J;% derivatives wrt reference coor
        Jac{j} = det(J);        if (Jac{j}<=0),error('Non-positive Jacobian');end
        dV=(Jac{j}*w(j));
        gradN_mean =gradN_mean+gradN{j}*dV;
        V =V+dV;
        self.phis(i)=max([self.phis(i),stab_fraction(self,J)]);
    end
    gradN_mean=gradN_mean/V; % Finish the calculation of the mean basis function gradient matrix
    %     Material orientation?
    if (~Rm_constant)% do I need to evaluate the local material orientation?
        c =mean(X);% physical location of the quadrature point
        if (~isempty(labels )),  Rm =Rmh(c,[],labels(i));%  No Jacobian matrix?
        else,                    Rm =Rmh(c,[],[]);                end
    end
    % Now we calculate the mean deformation gradient dx/dX   
    F1bar =x1'*gradN_mean;% wrt global material coordinates
    context.F =Rm'*F1bar*Rm;% Deformation gradient wrt  material orientation 
    % Update the stress for the real material
    [cauchy,self.matstates{i}] = update(self.material, self.matstates{i}, context);
    % Update the stress for the stabilization material
    %     [stabcauchy,~] = update(self.stabilization_material, [], context);
    %     Bbar = self.hBlmat(self,[],gradN_mean/F1bar,[],[]);% strain-displacement d/dx
    %     gcauchy =self.material.stress_vector_rotation(Rm')*(cauchy-f*stabcauchy); % to global
    %     Fe =  - Bbar'* (gcauchy * (V* det(F1bar))) ; % note the sign
    %     %     Now we update the stress for the stabilization material for the second time; this time for the full quadrature  rule
    %     for j=1:npts
    %         F1 =x1'*gradN{j};
    %         context.F = Rm'*F1*Rm;% Current deformation gradient wrt material orientation
    %         [stabcauchy,~] = update(self.stabilization_material, [], context);
    %         B = self.hBlmat(self,[],gradN{j}/F1,[],[]);% strain-displacement
    %         gcauchy =self.material.stress_vector_rotation(Rm')*(f*stabcauchy); % to global
    %         Fe = Fe - B'* (gcauchy * (Jac{j} * w(j)* det(F1))) ; % note the sign
    %     end
%     assemble(assembler, Fe, dofnums);
end
% F = make_vector (assembler);
return;
end