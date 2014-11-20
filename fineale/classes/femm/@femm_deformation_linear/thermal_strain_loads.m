function F = thermal_strain_loads(self, assembler, geom, u, dT)
% Compute the thermal-strain load vectors of the finite element set.
%
% function F = thermal_strain_loads(self, assembler, geom, u, dT)
%
% Return the assembled system vector F.
%    Arguments
%     assembler =  descendent of sysvec_assembler
%     geom=geometry field
%     u=displacement field
%     dT=temperature difference field (current temperature minus the
%         reference temperature at which the solid experiences no strains)
%
    fes = self.fes;% grab the finite elements to work on
    % Integration rule: compute the data needed for  numerical quadrature
    [npts Ns Nders w] = integration_data (self);
    % Material orientation
    Rm_constant = is_material_orientation_constant(self);% if not constant, need to compute  at each point
    if (~Rm_constant)
        Rmh = self.Rm;% handle to a function  to evaluate Rm
    else
        Rm = self.Rm;% constant material orientation matrix
    end    
    % Material
    mat = self.material;
    D_constant = are_tangent_moduli_constant (mat);
    % Retrieve data for efficiency
    conns = fes.conn; % connectivity
    labels = fes.label; % connectivity
    xs =geom.values;
    dTs =dT.values;
    % Prepare assembler
    Kedim =u.dim*fes.nfens;
    start_assembly(assembler, u.nfreedofs);
    % Now loop over all fes in the block
    for i=1:size(conns,1)
        conn =conns(i,:);
        x=xs(conn,:);
        DeltaT=dTs(conn,:);
        dofnums =reshape(u,gather_dofnums(u,conn));
        Fe =zeros(Kedim,1);;
        for j=1:npts
            c =Ns{j}'*x;% physical location of the quadrature point
            J = x' * Nders{j};% We compute the Jacobian matrix
            if (~Rm_constant)% do I need to evaluate the local material orientation?
                if (~isempty(labels )),  Rm =Rmh(c,J,labels(i));
                else,                    Rm =Rmh(c,J,[]);                end
            end
            Ndersp = Nders{j}/J;% derivatives wrt global coor
            Jac = Jacobian_volume(fes,conn, Ns{j}, J, x);
            if (Jac<=0),error('Non-positive Jacobian');end
            B = self.hBlmat(self,Ns{j},Ndersp*Rm,c,Rm);% strain-displacement
            context.dT =transpose(Ns{j})*DeltaT; context.xyz =c;
            iSigma = thermal_stress(mat,context);
            Fe = Fe - B'*(iSigma * (Jac * w(j)));
        end
        assemble(assembler, Fe, dofnums);
    end
    F = make_vector (assembler);
end

