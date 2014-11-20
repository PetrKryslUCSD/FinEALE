% Compute the conductivity matrix.
%
% function K = conductivity(self, assembler, geom, temp)
%
% Returns K as a matrix.
% Arguments
%           self  = heat diffusion model  
%           assembler = descendent of the sysmat_assembler class
%           geom=geometry field
%           temp=temperature field
function K = conductivity(self, assembler, geom, temp)
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
    start_assembly(assembler, Kedim, Kedim, size(conns,1), temp.nfreedofs, temp.nfreedofs);
    for i=1:size(conns,1)
        conn =conns(i,:);
        x=xs(conn,:);
        Ke =zeros(Kedim);;
        dofnums =reshape (temp,gather_dofnums(temp,conn));
        for j=1:npts
            J = Jacobian_matrix(fes,Nders{j},x);
            Jac = Jacobian_volume(fes,conn, Ns{j}, J, x);
            Ndersp = Nders{j}/J;
            if (~Rm_constant)% do I need to evaluate the local material orientation?
                if (~isempty(labels )),  Rm =Rmh(c,J,labels(i));
                else,                    Rm =Rmh(c,J,[]);                end
            end
            Ndersp=Ndersp*Rm;
            Ke = Ke + Ndersp*(kappa_bar*(Jac*w(j)))*Ndersp' ;
        end% Loop over quadrature points
        assembler.assemble_symmetric(Ke, dofnums);
    end% Loop over elements
    K = make_matrix (assembler);
end
