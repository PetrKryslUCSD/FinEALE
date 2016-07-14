function idat = inspect_integration_points(self, ...
        geom, un1, un, dt, dT, fe_list, context,...
        inspector, idat)
% Inspect the integration point quantities.
%
%     function idat = inspect_integration_points(self, ...
%             geom, u, dT, fe_list, context, inspector, idat)
%
% Input arguments
%    geom - reference geometry field
%    un1      - displacement field at the end of time step t_n+1; This is
%               the linear displacement field for linear problems
%    un       - displacement field at the end of time step t_n; This field
%               is ignored for linear problems
%    dt       - time step from  t_n to t_n+1; needed only by some
%                materials
%    dT - temperature difference field
%    fe_list - indexes of the finite elements that are to be inspected:
%          The fes to be included are: fes(fe_list).
%    context    - structure: see the update() method of the material.
%    inspector - function handle or in-line function with the signature
%             idat =inspector(idat, out, xyz, u, pc),
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
    if ~isfield(context,'output')
        context.output=  'cpress';% default output
    end
    % Material
    mat = self.material;
    ms=[];
    % Retrieve data for efficiency
    conns = fes.conn; % connectivity
    labels = fes.label; % finite element labels
    xs =geom.values;
    Us=un1.values;
    context.strain= [];
    if isempty(dT)
        dTs=zeros(geom.nfens,1);
    else
        dTs=dT.values;
    end
    % Now loop over selected fes 
    for m=1:length(fe_list)
        i=fe_list(m);
        conn = conns(i,:); % connectivity
        X=xs(conn,:);
        U=Us(conn,:);
        dT =dTs(conn,:);
        % Loop over all integration points
        for    qp =1:npts
            J = Jacobian_matrix(fes,Nders{qp},X);
            Jac = Jacobian_surface(fes,conn, Ns{qp}, J, X);
            qpx=Ns{qp}'*X;
            qpu=Ns{qp}'*U;
            [penetration,normal] = self.get_penetration(self.surface_data,qpx,qpu);
            cpress=0;
            if (penetration>0)
                cpress=(penetration)*self.penalty;%contact pressure
            end
            switch context.output
                case 'cpress'% contact pressure output
                    out =cpress;% To output coordinate system
            end
            if ~isempty (inspector)
                idat =feval(inspector,idat,out,Ns{qp}'*X,Ns{qp}'*U,pc(qp,:));
            end
        end; clear qp
    end
    return;
end
