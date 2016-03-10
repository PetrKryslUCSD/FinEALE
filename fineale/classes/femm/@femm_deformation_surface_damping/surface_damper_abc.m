% Compute the damper absorbing boundary condition matrices of the individual gcells.
% Return an array of them so they may be assembled.
%    Call as
%       C = surface_damper_abc(self, geom, v)
%    where
%       geom=geometry field
%       v=velocity field
%
function C = surface_damper_abc(self, assembler, geom, u)
   fes = self.fes;% grab the finite elements to work on
    % Integration rule: compute the data needed for  numerical quadrature
    [npts Ns Nders w] = integration_data (self);
    % Material orientation
    Rm_identity = is_material_orientation_identity(self);% if identity, work not needed
    Rm_constant = is_material_orientation_constant(self);% if not constant, need to compute  at each point
    if (~Rm_constant)
        Rmh = self.Rm;% handle to a function  to evaluate Rm
    else
        Rm = self.Rm;
    end    
    % Retrieve data for efficiency
    conns = fes.conn; % connectivity
    labels = fes.label; % connectivity
    xs =geom.values;
   % damping coefficient 
    h = self.damping_abc_impedance;
     % Prepare assembler
     dim=u.dim;
    Cedim =u.dim*fes.nfens;
    start_assembly(assembler, Cedim, Cedim, size(conns,1), u.nfreedofs, u.nfreedofs);
    % Now loop over all fes in the block
    for i=1:size(conns,1)
        conn =conns(i,:);
        x=xs(conn,:);
        dofnums =reshape(u,gather_dofnums(u,conn));
        Ce =zeros(Cedim);;
        for j=1:npts
            c =Ns{j}'*x;% physical location of the quadrature point
            J = x' * Nders{j};% We compute the Jacobian matrix
            if (~Rm_constant)% do I need to evaluate the local material orientation?
                if (~isempty(labels )),  Rm =Rmh(c,J,labels(i));
                else,                    Rm =Rmh(c,J,[]);                end
            end
            if (Rm_identity)
                Ndersp = Nders{j}/J;% derivatives wrt global coor
            else
                Ndersp = Nders{j}/(Rm'*J);% derivatives wrt local coor
            end
            Jac = Jacobian_surface(fes,conn, Ns{j}, J, x);
            if (Jac<=0),error('Non-positive Jacobian');end
            n=normal(self,c,J);
            for k= 1:fes.nfens 
                ns((k-1)*dim+1:k*dim)=n*Ns{j}(k);
            end 
            Ce = Ce + ns*ns' * (h* Jac * w(j));
        end
        assemble_symmetric(assembler, Ce, dofnums);
    end
    C = make_matrix (assembler);
end

