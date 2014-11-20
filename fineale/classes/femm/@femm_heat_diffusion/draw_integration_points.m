function draw_integration_points(self, gv, context)
% Produce graphic representation for material states in all finite elements.
%
% function draw_integration_points(self, gv, context)
%
%
% Input arguments
% self = FE model
% gv = graphic viewer
% context = struct
% with typically mandatory fields
%    x=reference geometry field
%    theta=temperature field
% with optional fields
%    quantity = (default 'flux')
%    scale = scale the representation of the material state graphic by this
%       factor
% Any other fields in the context struct may be interpreted by the
% implementation of the material: see the draw() function for the
% particular material.
%
    quantity='flux';
    if (isfield(context,'quantity'))
        quantity= context.quantity;
    end
    scale = 1.0;
    if (isfield(context,'scale'))
        scale= context.scale;
    end
    function idat=inspector(idat, out, xyz, u, pc)
        switch (quantity)
            case 'flux'
                flux = out;
                component =1;
                if isfield(context,'component')
                    component = context.component;
                end
                color =[1 1 0];
                if isfield(context,'data_cmap')
                    color = map_data(context.data_cmap, flux(component));
                end
                context.facecolor= color;
                draw_arrow(gv, xyz, scale*flux, context);
            otherwise
                %  do nothing;
        end
    end
    idat= [];
    idat = inspect_integration_points(self, ...
                    context.x, context.theta, 1:count(self.fes), context,...
                    @inspector, idat);
end
    % fes = self.fes;
    % dim = context.x.dim;
    % % Precompute basis f. values + basis f. gradients wrt parametric coor
    % [npts Ns Nders w] = integration_data (self);
    % % Material
    % mat = self.material;
    % matstates=self.matstates;
    % % Note that the thermal conductivity matrix is in the 
    % % local  material orientation coordinates.
    % kappa_bar =  mat.property.thermal_conductivity;
    % % Material orientation?
    % Rm_identity = is_material_orientation_identity(self);
    % Rm_constant = is_material_orientation_constant(self);
    % if (~Rm_identity) && (Rm_constant)
        % Rm = self.Rm;
    % end   
    % % Prepare some data: 
    % conns = fes.conn; % connectivity
    % xs =context.x.values;% retrieve the geometry information
    % Ts =context.theta.values;% retrieve the geometry information
    % context.update_context= [];
    % for i=1:size(conns,1)
        % conn =conns(i,:);
        % x=xs(conn,:);
        % T=Ts(conn,:);
        % for j=1:npts
            % J = Jacobian_matrix(fes,Nders{j},x);
            % Jac = Jacobian_volume(fes,conn, Ns{j}, J, x);
            % Ndersp = Nders{j}/J;
            % if (~Rm_identity)
                % if (~Rm_constant)
                    % Rm = material_orientation(self,fes,pc(j,:),x);
                % end
                % Ndersp=Ndersp*Rm;% gradients to material orientation coordinates
                % context.Rm=Rm;% the draw method needs to  transform
            % end  
            % context.ms=matstates{i,j};
            % context.xyz=transpose(Ns{j})*x;
            % context.update_context.gradtheta = transpose(T)* Ndersp;
            % draw(mat, gv, context);
        % end% Loop over quadrature points
    % end% Loop over elements
% end
