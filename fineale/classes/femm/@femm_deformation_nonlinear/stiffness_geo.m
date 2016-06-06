function K = stiffness_geo (self, assembler, geom, un1, un, dt)
% Compute the geometric stiffness matrix by computing and assembling the
% matrices of the individual FEs.
%
% function K = stiffness_geo(self, assembler, geom, u1, u)
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
% Integration rule
[npts Ns gradNs w] = integration_data (self);;
% Material orientation
Rm_constant = is_material_orientation_constant(self);% if not constant, need to compute  at each point
if (~Rm_constant)
    Rmh = self.Rm;% handle to a function  to evaluate Rm
else
    Rm = self.Rm;
end
%  Indexing vector
idx = (1:un1.dim:(self.fes.nfens-1)*un1.dim+1);
% Retrieve data for efficiency
conns = self.fes.conn; % connectivity
labels = self.fes.label; % connectivity
Xs = geom.values; % reference coordinates
Uns = un.values; % displacement in step n
Un1s = un1.values; % displacement in step n+1
context.Fn1 = []; context.Fn = []; context.dt=dt;
context.dT = [];
context.output='Cauchy';
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
    K1 =zeros(Kedim);;
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
        Fn1 = xn1'*gradNX;% Current deformation gradient
        context.Fn1=Rm'*Fn1*Rm;%  deformation gradient  in material coordinates
        Fn = xn'*gradNX;% Current deformation gradient
        context.Fn=Rm'*Fn*Rm;%  deformation gradient  in material coordinates
        [cauchy,self.matstates{i,j}] = update(self.material, self.matstates{i,j}, context);
        gcauchy =self.material.stress_vector_rotation(Rm')*cauchy; % material to global
        sigma = stress_6v_to_3x3t(self.material,gcauchy);
        gradNx1 =gradNX/Fn1;
        c1 = gradNx1*sigma*gradNx1';
        K1(idx,idx)     = c1;
        K1(idx+1,idx+1) = c1;
        K1(idx+2,idx+2) = c1;
        Ke = Ke + K1*(Jac*w(j)*det(Fn1));
    end
    assemble_symmetric(assembler, Ke, dofnums);
end
K = make_matrix (assembler);
end
