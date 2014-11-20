function elerrs = flux_L2_error (self, geom, u, dT, nodal_stress)
% Compute the L2 stress flux error of the individual gcells.
%
% function elerrs = flux_L2_error (self, geom, u, dT, nodal_stress)
%
% Input arguments
%     self     - finite element block
%     geom     - reference geometry field
%     u        - displacement field
%     dT       - temperature difference field, may be supplied as empty
%     nodal_stress   -  field of the nodal values of the stress. This would be computed as
%       nodal_stress = field_from_integration_points_spr(femm, geom, u, [], 'Cauchy', 1:4);
%       for a plane strain problem (as an example).  
% 
% Output argument
% Return an array of the element L2 errors of the flux.
%     elerrs= Array, number of columns equals the number of geometric
%     cells; for each finite element e the array component elerrs(e) holds the 
%     L2 norm of the error of the flux over the element 
%      
    fes = self.fes;% grab the finite elements to work on
    % Integration rule: compute the data needed for  numerical quadrature
    [npts Ns Nders w pc] = integration_data (self);
    % Material orientation matrix
    Rm_constant = is_material_orientation_constant(self);% if not constant, need to compute  at each point
    if (~Rm_constant)
        Rmh = self.Rm;% handle to a function  to evaluate Rm
    else
        Rm = self.Rm;% constant material orientation matrix
    end    
    % Material
    mat = self.material;
    % Retrieve data for efficiency
    conns = fes.conn; % connectivity
    labels = fes.label; % finite element labels
    xs =geom.values;
    if isempty(dT)
        dTs=zeros(geom.nfens,1);
    else
        dTs=dT.values;
    end
    % output
    dim=nodal_stress.dim;
    errors=zeros(dim,count(fes));
    % Now loop over all gcells in the block
    for i=1:count(fes)
        conn = conns(i,:); % connectivity
        x=xs(conn,:);
        U=reshape(u, gather_values(u,conn));
        dT =dTs(conn,:);
        nodsig = gather_values(nodal_stress, conn);
         % Loop over all integration points
        for j=1:npts
            c =Ns{j}'*x;% physical location
            J = x' * Nders{j};
            if (~Rm_constant)% do I need to evaluate the local material orientation?
                if (~isempty(labels )),  Rm =Rmh(c,J,labels(i));
                else,                    Rm =Rmh(c,J,[]);                end
            end
            Ndersp = Nders{j}/J;% derivatives wrt global coor
            Jac = Jacobian_volume(fes,conn, Ns{j}, J, x);
            B = self.hBlmat(self,Ns{j},Ndersp*Rm,c,Rm);% strain-displacement
            context.strain = B*U;%strain in the local material orientation coord
            context.dT = transpose(Ns{j})*dT;
            context.xyz =c;
            [sig,ignore] = update(mat, [], context);
            errors(:,i)=errors(:,i)+(sig-(Ns{j}'*nodsig)').^2 * Jac * w(j);
        end
    end
    elerrs=sqrt(sum(errors,1));
end
