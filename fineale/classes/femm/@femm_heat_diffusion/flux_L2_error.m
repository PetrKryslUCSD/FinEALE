% Compute the  L2 heat flux error for the individual gcells.
%
% function elerrs = flux_L2_error (self, geom, Temp, nodal_flux)
%
% Input arguments
%     self     - finite element block
%     geom     - reference geometry field
%     u        - displacement field
%     dT       - temperature difference field, may be supplied as empty
%     nodal_flux   -  field of the nodal values of the heat flux. This would be computed as
%       nodal_flux = field_from_integration_points(femm, geom, u, [], 'flux', 1:2);
%       for a plane strain problem (as an example).  
% 
% Output argument
% Return an array of the element L2 errors of the flux.
%     elerrs= Array, number of columns equals the number of geometric
%     cells; for each finite element e the array component elerrs(e) holds the 
%     L2 norm of the error of the flux over the element 
%      
function elerrs = flux_L2_error (self, geom, Temp, nodal_flux)
    fes = self.fes;% grab the finite elements to work on
    % Integration rule: compute the data needed for  numerical quadrature
    [npts Ns Nders w pc] = integration_data (self);
    % Material
    mat = self.material;
    matstates=self.matstates;
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
    % output
    dim=nodal_flux.dim;
    errors=zeros(dim,count(fes));
    conns = fes.conn; % connectivity
    xs =geom.values;% retrieve the geometry information
    Ts =Temp.values;% retrieve the geometry information
    nodsigs =nodal_flux.values;% retrieve the geometry information
    % Now loop over all gcells in the block
    for i=1:size(conns,1)
        conn = conns(i,:); % connectivity
        x = xs(conn, :); % coord
        T = Ts(conn); % temperature
        nodsig = nodsigs(conn, :); % displ
         % Loop over all integration points
        for j=1:npts
            J = Jacobian_matrix(fes,Nders{j},x);
            Jac = Jacobian_volume(fes,conn, Ns{j}, J, x);
            if (~Rm_constant)% do I need to evaluate the local material orientation?
                if (~isempty(labels )),  Rm =Rmh(c,J,labels(i));
                else,                    Rm =Rmh(c,J,[]);                end
            end
            Ndersp = Nders{j}/(Rm'*J);% gradient WRT the material coordinates
            context.xyz=Ns{j}'*x;
            context.gradtheta = T'* Ndersp;
            [sig,ignore] = update(mat, matstates{i,j}, context);
            errors(:,i)=errors(:,i)+(sig-(Ns{j}'*nodsig)').^2 * Jac * w(j);
        end
    end
    elerrs=sqrt(sum(errors,1));
end
