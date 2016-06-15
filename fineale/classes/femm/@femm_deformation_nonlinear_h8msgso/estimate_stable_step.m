function stabldt = estimate_stable_step (self, geom)
% Estimate the stable time step for explicit time integration.  
%
% function stabldt = estimate_stable_step (self, geom, un1, un, dt)
%
% Return an estimate of the stable time step.
%     geom=geometry field
%     un1      - displacement field at the end of time step t_n+1
%     un       - displacement field at the end of time step t_n
%     dt       - time step from  t_n to t_n+1; needed only by some
%                materials

if isempty(self.phis)%  Do we need to calculate the form factors?
    error('Need the stabilization parameters');
end

% Integration rule
[npts, Ns, gradNparams, w] = integration_data (self);;
gradN =cell(8,1); Jac=cell(8,1);
% Retrieve data for efficiency
conns = self.fes.conn; % connectivity
Xs = geom.values; % reference coordinates
dsq=@(X1,X2)sum((X2-X1).*(X2-X1));
stabldt=inf;
for i=1:size(conns,1)
    conn =conns(i,:); % connectivity
    X=Xs(conn,:);
    mind =sqrt(min([...
        dsq(X(1,:),X(2,:)),...
        dsq(X(2,:),X(3,:)),...
        dsq(X(3,:),X(4,:)),...
        dsq(X(4,:),X(1,:)),...
        dsq(X(1,:),X(5,:)),...
        dsq(X(2,:),X(6,:)),...
        dsq(X(3,:),X(7,:)),...
        dsq(X(4,:),X(8,:)),...
        dsq(X(5,:),X(6,:)),...
        dsq(X(6,:),X(7,:)),...
        dsq(X(7,:),X(8,:)),...
        dsq(X(8,:),X(5,:))]));
    speed_of_sound = sqrt(self.material.property.E / (1 - 2*self.material.property.nu) / self.material.property.rho);
    stabldt = min([stabldt, mind/speed_of_sound]);
end


end
