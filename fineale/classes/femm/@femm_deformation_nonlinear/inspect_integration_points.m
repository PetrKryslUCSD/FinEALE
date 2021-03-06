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
    % Material orientation matrix
    Rm_constant = is_material_orientation_constant(self);% if not constant, need to compute  at each point
    if (~Rm_constant)
        Rmh = self.Rm;% handle to a function  to evaluate Rm
    else
        Rm = self.Rm;
    end    
    % Output orientation matrix
    outputRm_identity = ~isfield(context,'outputRm');% if identity, work not needed
    outputRm_constant = (outputRm_identity) ...
        || strcmp(class(context.outputRm),'double');% need to compute at each point
    if (~outputRm_constant)
        outputRmh = context.outputRm;% handle to a function  to evaluate outputRm
    else
        if (~outputRm_identity)
            outputRm = context.outputRm;
        end
    end    
    if ~isfield(context,'output')
        context.output=  'Cauchy';% default output
    end
    % Material
    mat = self.material;
    % Retrieve data for efficiency
    conns = fes.conn; % connectivity
    labels = fes.label; % finite element labels
    Xs =geom.values;
    Uns = un.values; % displacement in step n
    Un1s = un1.values; % displacement in step n+1
    context.Fn1 = []; context.Fn = []; context.dt=dt;
    if isempty(dT)
        dTs=zeros(geom.nfens,1);
    else
        dTs=dT.values;
    end
    % Now loop over selected fes 
    for m=1:length(fe_list)
        i=fe_list(m);
        conn = conns(i,:); % connectivity
        conn =conns(i,:);
        X=Xs(conn,:);
        Un=Uns(conn,:);
        Un1=Un1s(conn,:);
        xn = X + Un; % previous known coordinates
        xn1 = X + Un1; % current coordinates
        dT =dTs(conn,:);
        % Loop over all integration points
        for j=1:npts
            c =Ns{j}'*X;% physical location
            u_c = transpose(Ns{j})*Un1s(conn,:);% Current displacement of quadrature point
            J = X' * Nders{j};
            if (~Rm_constant)% do I need to evaluate the local material orientation?
                if (~isempty(labels )),  Rm =Rmh(c,J,labels(i));
                else,                    Rm =Rmh(c,J,[]);                end
            end
            if (~outputRm_constant)% do I need to evaluate the output orientation?
                if (~isempty(labels )),  outputRm =outputRmh(c,J,labels(i));
                else,                    outputRm =outputRmh(c,J,[]);    end
            end
            gradNX = Nders{j}/J;% derivatives wrt global coor
            Jac = Jacobian_volume(fes,conn, Ns{j}, J, X);
            Fn1 = xn1'*gradNX;% Current deformation gradient
            context.Fn1=Rm'*Fn1*Rm;%  deformation gradient  in material coordinates
            Fn = xn'*gradNX;% Current deformation gradient
            context.Fn=Rm'*Fn*Rm;%  deformation gradient  in material coordinates
            context.dT = transpose(Ns{j})*dT;
            context.xyz =c;
            context.ms=self.matstates{i,j};
            out = state(mat, context.ms, context);
            switch context.output
                case 'Cauchy'
                    if (~outputRm_constant)
                        Rmc =mean(X);% physical location of the quadrature point
                        if (~isempty(labels )),  outputRm=outputRmh(Rmc,[],labels(i));%  No Jacobian matrix?
                        else,                    outputRm =outputRmh(Rmc,[],[]);                end
                    end
                    if (outputRm_identity)
                        out =mat.stress_vector_rotation((Rm'))*out;% To output coordinate system
                    else % From material, to global, to output
                        out =mat.stress_vector_rotation((Rm'*outputRm))*out;% To output coordinate system
                    end
            end
            if ~isempty (inspector)
                idat =feval(inspector,idat,out,c,u_c,pc(j,:));
            end
        end
    end
    return;
end

