function K = stiffness_geo (self, assembler, geom, u1,u)
% Compute the geometric stiffness matrix by computing and assembling the
% matrices of the individual FEs.
%
% function K = stiffness_geo(self, assembler, geom, u1, u)
%
% Return a stiffness matrix.
%     assembler = descendent of the sysmat_assembler class
%     geom=geometry field
%     u1=displacement field, the current configuration
%     u=displacement field, the last known configuration

% Integration rule
[npts, Ns, gradNparams, w] = integration_data (self);;
gradN =cell(8,1); Jac=cell(8,1);
% Material orientation
Rm_constant = is_material_orientation_constant(self);% if not constant, need to compute  at each point
if (~Rm_constant)
    Rmh = self.Rm;% handle to a function  to evaluate Rm
else
    Rm = self.Rm;
end
%  Indexing vector
idx = (1:u.dim:(self.fes.nfens-1)*u.dim+1);
% Retrieve data for efficiency
conns = self.fes.conn; % connectivity
labels = self.fes.label; % connectivity
Xs = geom.values; % reference coordinates
Us = u.values; % displacement in step n
U1s = u1.values; % displacement in step n+1
context.dT = [];
context.F = [];
context.output='Cauchy';
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
    K1 =zeros(Kedim);;
    % First we calculate  the mean basis function gradient matrix and the volume of the element
    gradN_mean =zeros(self.fes.nfens,u.dim); V=0;
    for j=1:npts % loop over all quadrature points
        J = X' * gradNparams{j};% Jacobian matrix wrt reference coordinates
        gradN{j} = gradNparams{j}/J;% derivatives wrt global coor
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
    F1bar =x1' * gradN_mean;% With respect to reference coordinates
    gradN_x1_mean= gradN_mean/F1bar; %...  and the b.f. gradients with respect to current conf
    context.F =Rm'*F1bar*Rm;% Deformation gradient wrt  material orientation 
    [cauchy,~] = update(self.material, self.matstates{i}, context);
    [stabcauchy,~] = update(self.stabilization_material, [], context);
    gcauchy =self.material.stress_vector_rotation(Rm')*(cauchy-self.phis(i)*stabcauchy); % to global
    sigma = stress_6v_to_3x3t(self.material,gcauchy);
    c1 = gradN_x1_mean*sigma*gradN_x1_mean';
    K1(idx,idx)     = c1;
    K1(idx+1,idx+1) = c1;
    K1(idx+2,idx+2) = c1;
    Ke = K1*V;
    for j=1:npts       % Loop over all integration points
        F1 = x1'*gradN{j};% Current deformation gradient
        context.F = Rm'*F1*Rm;% Current deformation gradient wrt material orientation
        [stabcauchy,~] = update(self.stabilization_material, [], context);
        gcauchy =self.material.stress_vector_rotation(Rm')*(self.phis(i)*stabcauchy); % to global
        sigma = stress_6v_to_3x3t(self.stabilization_material,gcauchy);
        gradNx1 =gradN{j}/F1;
        c1 = gradNx1*sigma*gradNx1';
        K1(idx,idx)     = c1;
        K1(idx+1,idx+1) = c1;
        K1(idx+2,idx+2) = c1;
        Ke = Ke + K1*(Jac{j}*w(j));
    end
    assemble_symmetric(assembler, Ke, dofnums);
end
K = make_matrix (assembler);
end
