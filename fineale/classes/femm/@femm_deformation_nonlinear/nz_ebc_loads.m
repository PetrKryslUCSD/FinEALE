function F = nz_ebc_loads(self, assembler, geom, un1, un, du, dt)
% Compute load vector for nonzero essential boundary conditions.
%
% function F = nz_ebc_loads(self, assembler, geom, un1, un, dt)
%
% Return the load vector matrix corresponding to prescribed displacements.
%     assembler = descendent of the sysvec_assembler class
%     geom=geometry field
%     un1      - displacement field at the end of time step t_n+1
%     un       - displacement field at the end of time step t_n
%     du       - displacement increment imposed at t_n+1
%     dt       - time step from  t_n to t_n+1; needed only by some
%                materials
%

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
start_assembly(assembler, un1.nfreedofs);
for i=1:size(conns,1)
    conn =conns(i,:); % connectivity
    pu = reshape(du, gather_fixed_values(du, conn));% Retrieve the displacement increment
    if norm (pu) ~= 0
        X=Xs(conn,:);
        Un=Uns(conn,:);
        Un1=Un1s(conn,:);
        xn = X + Un; % previous known coordinates
        xn1 = X + Un1; % current coordinates
        dofnums =reshape(un1,gather_dofnums(un1,conn));
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
            Fn1 = xn1'*gradNX;% Current deformation gradient
            context.Fn1=Rm'*Fn1*Rm;%  deformation gradient  in material coordinates
            Fn = xn'*gradNX;% Current deformation gradient
            context.Fn=Rm'*Fn*Rm;%  deformation gradient  in material coordinates
            B = self.hBlmat(self,Ns{j},gradNX/Fn1,c,[]);% strain-displacement
            context.ms =self.matstates{i,j};
            D = tangent_moduli (self.material, context);% in material coordinates
            D = rotate_stiffness(self.material,D,Rm');%in global coordinates
            Ke = Ke + (B'*(D*(Jac*w(j)*det(Fn1)))*B);
        end
        Fe =  -Ke*pu;
        assemble(assembler, Fe, dofnums);
    end
end
F = make_vector (assembler);
end
