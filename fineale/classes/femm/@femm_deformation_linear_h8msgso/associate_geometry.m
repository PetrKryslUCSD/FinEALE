function self=associate_geometry(self,geom)
% Associate a geometry field with the FEMM.
%
%         function self=associate_geometry(self,geom)
%
% geom= Geometry field to associate with self. Pass in []
%    (empty array) to disassociate any existing geometry from the
%    self FEMM.
%
% Here we use the association with the geometry to compute the
% stabilization factors based on the material properties and the shapes of
% the elements.
associate_geometry@femm_base(self,geom);
% Now compute and store the stabilization parameters
[npts, Ns, gradNparams, w] = integration_data (self);;
% Retrieve data for efficiency
conns = self.fes.conn; % connectivity
labels = self.fes.label; % connectivity
Xs = geom.values; % reference coordinates
self.phis=zeros(size(conns,1),1);
for i=1:size(conns,1)
    conn =conns(i,:); % connectivity
    X=Xs(conn,:);
    % First we calculate  the mean basis function gradient matrix and the volume of the element
    V=0;
    self.phis(i)=0;%  stabilization fraction
    for j=1:npts % loop over all quadrature points
        J = X' * gradNparams{j};% Jacobian matrix wrt reference coordinates
        Jac{j} = det(J);        if (Jac{j}<=0),error('Non-positive Jacobian');end
        dV=(Jac{j}*w(j));
        V =V+dV;
        self.phis(i)=max([self.phis(i),stab_fraction(self,J)]);
    end
end
end
