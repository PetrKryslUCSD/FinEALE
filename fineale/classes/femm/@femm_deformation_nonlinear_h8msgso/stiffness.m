
function K = stiffness (self, assembler, geom, un1, un, dt)
% Compute the stiffness matrix by computing and assembling the
% matrices of the individual FEs.
%
% function K = stiffness (self, assembler, geom, un1, un, dt)
%
% Return a stiffness matrix.
%     assembler = descendent of the sysmat_assembler class
%     geom=geometry field
%     un1      - displacement field at the end of time step t_n+1
%     un       - displacement field at the end of time step t_n
%     dt       - time step from  t_n to t_n+1; needed only by some
%                materials

if (~exist('dt','var'))
    % If the time step is needed, this will cause trouble.  Which indicates
    % it should have been supplied. 
    dt=[];
end
if (~exist('un','var'))
    % When this method gets called in linear problems, no previous
    % deformation is supplied: make one up.
    un=0*un1;
end
if isempty(self.phis)%  Do we need to calculate the form factors?
    error('Need the stabilization parameters');
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
Uns = un.values; % displacement in step n
Un1s = un1.values; % displacement in step n+1
context.stiff_type='Eulerian';
context.dt=dt;% time step
% Prepare assembler
Kedim =un1.dim*self.fes.nfens;
start_assembly(assembler, Kedim, Kedim, size(conns,1), un1.nfreedofs, un1.nfreedofs);
for i=1:size(conns,1)
    conn =conns(i,:); % connectivity
    X=Xs(conn,:);
    Un=Uns(conn,:);
    Un1=Un1s(conn,:);
    xn = X + Un; % previous known coordinates
    xn1 = X + Un1; % current coordinates
    dofnums =reshape(un1,gather_dofnums(un1,conn));
    Ke =zeros(Kedim);;
    % First we calculate  the mean basis function gradient matrix and the volume of the element
    gradN_mean =zeros(self.fes.nfens,un1.dim); V=0;
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
    Fnbar =xn' * gradN_mean;
    Fn1bar =xn1' * gradN_mean;
    Bbar = self.hBlmat(self,[],gradN_mean/Fn1bar,[],[]);% strain-displacement d/dX
    context.ms=self.matstates{i};
    context.Fn=Rm'*Fnbar*Rm;%  Deformation gradient in material coordinates
    context.Fn1=Rm'*Fn1bar*Rm;%  Deformation gradient in material coordinates
    D = tangent_moduli (self.material, context);
    D = rotate_stiffness(self.material,D,Rm');
    Ke = (Bbar'*(D*(V))*Bbar);
    context.ms=[]; % No mat. state for stabil. material needed: only Hyperelastic
    Dstab = tangent_moduli (self.stabilization_material, context);
    Dstab = rotate_stiffness(self.stabilization_material,Dstab,Rm');
    Ke = Ke - (Bbar'*(Dstab*(self.phis(i)*V))*Bbar);
    for j=1:npts       % Loop over all integration points
        Fn1 = xn1'*gradN{j};% Current deformation gradient
        B = self.hBlmat(self,[],gradN{j}/Fn1,[],[]);% strain-displacement
        context.Fn1=Rm'*Fn1*Rm;
        Dstab = tangent_moduli (self.stabilization_material, context);
        Dstab = rotate_stiffness(self.stabilization_material,Dstab,Rm');
        Ke = Ke + (B'*(Dstab*(self.phis(i)*Jac{j}*w(j)))*B);
    end
    assemble_symmetric(assembler, Ke, dofnums);
end
K = make_matrix (assembler);
end


function K = stiffness_for_anisotropic_stabilization_material (self, assembler, geom, un1, un, dt)
% Compute the stiffness matrix by computing and assembling the
% matrices of the individual FEs.
%
% function K = stiffness (self, assembler, geom, un1, un, dt)
%
% Return a stiffness matrix.
%     assembler = descendent of the sysmat_assembler class
%     geom=geometry field
%     un1      - displacement field at the end of time step t_n+1
%     un       - displacement field at the end of time step t_n
%     dt       - time step from  t_n to t_n+1; needed only by some
%                materials

if (~exist('dt','var'))
    % If the time step is needed, this will cause trouble.  Which indicates
    % it should have been supplied. 
    dt=[];
end
if (~exist('un','var'))
    % When this method gets called in linear problems, no previous
    % deformation is supplied: make one up.
    un=0*un1;
end
if isempty(self.phis)%  Do we need to calculate the form factors?
    error('Need the stabilization parameters');
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
Uns = un.values; % displacement in step n
Un1s = un1.values; % displacement in step n+1
context.stiff_type='Eulerian';
context.dt=dt;% time step
% Prepare assembler
Kedim =un1.dim*self.fes.nfens;
start_assembly(assembler, Kedim, Kedim, size(conns,1), un1.nfreedofs, un1.nfreedofs);
for i=1:size(conns,1)
    conn =conns(i,:); % connectivity
    X=Xs(conn,:);
    Un=Uns(conn,:);
    Un1=Un1s(conn,:);
    xn = X + Un; % previous known coordinates
    xn1 = X + Un1; % current coordinates
    dofnums =reshape(un1,gather_dofnums(un1,conn));
    Ke =zeros(Kedim);;
    % First we calculate  the mean basis function gradient matrix and the volume of the element
    gradN_mean =zeros(self.fes.nfens,un1.dim); V=0;
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
    Fnbar =xn' * gradN_mean;
    Fn1bar =xn1' * gradN_mean;
    Bbar = self.hBlmat(self,[],gradN_mean/Fn1bar,[],[]);% strain-displacement d/dX
    context.ms=self.matstates{i};
    context.Fn=Rm'*Fnbar*Rm;%  Deformation gradient in material coordinates
    context.Fn1=Rm'*Fn1bar*Rm;%  Deformation gradient in material coordinates
    D = tangent_moduli (self.material, context);
    D = rotate_stiffness(self.material,D,Rm');
    Ke = (Bbar'*(D*(V))*Bbar);
    context.ms=[]; % No mat. state for stabil. material needed: only Hyperelastic
    Dstab = tangent_moduli (self.stabilization_material, context);
    Dstab = rotate_stiffness(self.stabilization_material,Dstab,Rm');
    Ke = Ke - (Bbar'*(Dstab*(self.phis(i)*V))*Bbar);
    for j=1:npts       % Loop over all integration points
        Fn1 = xn1'*gradN{j};% Current deformation gradient
        B = self.hBlmat(self,[],gradN{j}/Fn1,[],[]);% strain-displacement
        context.Fn1=Fn1;
        Dstab = tangent_moduli (self.stabilization_material, context);
        Ke = Ke + (B'*(Dstab*(self.phis(i)*Jac{j}*w(j)))*B);
    end
    assemble_symmetric(assembler, Ke, dofnums);
end
K = make_matrix (assembler);
end
