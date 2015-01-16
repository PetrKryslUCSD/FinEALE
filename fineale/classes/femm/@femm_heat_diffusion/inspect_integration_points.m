function idat = inspect_integration_points(self, ...
        geom, Temp, gcell_list, context,...
        inspector, idat)
% Inspect the integration point quantities.
%
% function idat = inspect_integration_points(self, ...
%         geom,  Temp, gcell_list, context,...
%         inspector, idat)
%
% Input arguments
%    geom - reference geometry field
%    Temp - temperature field
%    gcell_list - indexes of the geometric cells that are to be inspected:
%          The gcells to be included are: gcells(gcell_list).
%    context    - structure: see the update() method.
%    inspector - function handle or in-line function with the signature
%             idat =inspector(idat, out, xyz, pc),
%        where
%         idat - a structure or an array that the inspector may
%                use to maintain some state, for instance minimum or
%                maximum of stress, out is the output  of the update()
%                method, xyz is the location of the integration point
%                in the *reference* configuration, and u the displacement 
%                of the integration point. The argument pc are the
%                parametric coordinates of the quadrature point.
% Output arguments
%     idat - see the description of the input argument
%
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
    % Prepare some data: 
    conns = fes.conn; % connectivity
    xs =geom.values;% retrieve the geometry information
    Ts =Temp.values;% retrieve the geometry information
   % Now loop over all gcells in the block
    for m=1:length(gcell_list)
        i=gcell_list(m);
        conn =conns(i,:);
        x = xs(conn, :); % coord
        T = Ts(conn); % temperature
        % Loop over all integration points
        for j=1:npts
            J = Jacobian_matrix(fes,Nders{j},x);
            Jac = Jacobian_volume(fes,conn, Ns{j}, J, x);
            if (~Rm_constant)% do I need to evaluate the local material orientation?
                if (~isempty(labels )),  Rm =Rmh(c,J,labels(i));
                else,                    Rm =Rmh(c,J,[]);                end
            end
            Ndersp = Nders{j}/(Rm'*J);% gradient WRT the material coordinates
            context.xyz=Ns{j}'*x;% 'location in global coordinates
            context.gradtheta = T'* Ndersp;%  gradient in mat. orient. coord
            context.Rm  =Rm;% material orientation matrix
            % Now use the material update method to compute the requested
            % quantity ...
            [out,ignore] = update(mat, matstates{i,j}, context);
            %  ... and call the callback to allow for its processing
            if ~isempty (inspector)
                idat =feval(inspector,idat,out,context.xyz,[],pc(j,:));
            end
        end
    end
    return;
end
