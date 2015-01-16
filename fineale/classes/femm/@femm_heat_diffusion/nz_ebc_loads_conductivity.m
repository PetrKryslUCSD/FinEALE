function F = nz_ebc_loads_conductivity(self, assembler, geom, temp)
% Compute load vector for nonzero EBC for fixed temperature.
%
% function F = nz_ebc_loads_conductivity(self, assembler, geom, temp)
%
%    Arguments
%      assembler =  descendent of sysvec_assembler
%      geom=geometry field
%      temp=temperature field
%
% Return the assembled system vector F.

fes = self.fes;
    % Precompute basis f. values + basis f. gradients wrt parametric coor
    [npts Ns Nders w] = integration_data (self);
    % Material
    mat = self.material;
    % Note that the thermal conductivity matrix is in the 
    % local  material orientation coordinates.
    kappa_bar =  mat.property.thermal_conductivity;
    % Material orientation?
    Rm_constant = is_material_orientation_constant(self);
    if (~Rm_constant)
        Rmh = self.Rm;% handle to a function  to evaluate Rm
    else
        Rm = self.Rm;% constant material orientation matrix
    end    
    % Prepare some data: 
    conns = fes.conn; % connectivity
    labels = fes.label; % connectivity
    xs =geom.values;% retrieve the geometry information
    % Prepare assembler
    Kedim =temp.dim*fes.nfens;
    start_assembly(assembler, temp.nfreedofs);
    for i=1:size(conns,1)
        conn =conns(i,:);
        pT = reshape(temp, gather_fixed_values(temp, conn));
        if norm (pT) ~= 0
            x=xs(conn,:);
            Ke =zeros(Kedim);;
            dofnums =reshape (temp,gather_dofnums(temp,conn));
            for j=1:npts
                J = Jacobian_matrix(fes,Nders{j},x);
                Jac = Jacobian_volume(fes,conn, Ns{j}, J, x);
                if (~Rm_constant)% do I need to evaluate the local material orientation?
                    if (~isempty(labels )),  Rm =Rmh(c,J,labels(i));
                    else,                    Rm =Rmh(c,J,[]);                end
                end
                Ndersp = Nders{j}/(Rm'*J);% gradient WRT the material coordinates
                Ke = Ke + Ndersp*(kappa_bar*(Jac*w(j)))*Ndersp' ;
            end% Loop over quadrature points
            Fe =  -Ke*pT;
            assemble(assembler, Fe, dofnums);
        end
    end% Loop over elements
    F= make_vector (assembler);
end
