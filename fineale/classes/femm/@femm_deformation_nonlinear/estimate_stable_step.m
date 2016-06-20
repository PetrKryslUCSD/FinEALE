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


% Retrieve data for efficiency
conns = self.fes.conn; % connectivity
Xs = geom.values; % reference coordinates
dsq=@(X1,X2)sum((X2-X1).*(X2-X1));
stabldt=inf;
for i=1:size(conns,1)
    conn =conns(i,:); % connectivity
    X=Xs(conn,:);
    mind =inf;
    for  i =1:length(conn)
        for  j =(i+1):length(conn)
            mind=min([mind,dsq(X(i,:),X(j,:))]);
        end
    end
    mind =sqrt(mind);
    % The formula below is strictly speaking only applicable in three
    % dimensions for isotropic materials.  This needs to be generalized at
    % some point.
    speed_of_sound = sqrt(self.material.property.E / (1 - 2*self.material.property.nu) / self.material.property.rho);
    stabldt = min([stabldt, mind/speed_of_sound]);
end


end
