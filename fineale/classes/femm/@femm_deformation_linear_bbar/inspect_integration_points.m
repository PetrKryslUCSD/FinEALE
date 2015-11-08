function idat = inspect_integration_points(self, ...
        geom, u, dT, fe_list, context,...
        inspector, idat)
% Inspect the integration point quantities.
%
% function idat = inspect_integration_points(self, ...
        % geom, u, dT, fe_list, context,...
        % inspector, idat)
%
% Input arguments
%    geom - reference geometry field
%    u - displacement field
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
    [npts_c,Ns_c,Nders_c,Ns_pv_c,w_c,pc_c] = integration_data_constrained (self);
    [npts_u,Ns_u,Nders_u,Ns_pv_u,w_u,pc_u] = integration_data_unconstrained (self);
    % Material orientation
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
    D_constant = are_tangent_moduli_constant (mat);
    if (D_constant)
        D = tangent_moduli(mat,[]);
        [m1,Id,Iv]=constrained_mat_data(self,D);
    end
    PVdim=size(Ns_pv_c{1},1);
    Ddim=size(D,1);;
    Kedim =u.dim*fes.nfens;
    % Retrieve data for efficiency
    conns = fes.conn; % connectivity
    xs =geom.values;
    Us =u.values;
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
        x=xs(conn,:);
        U=reshape(u, gather_values(u,conn));
        dT =dTs(conn,:);
        % Loop over all integration points
        CT =zeros(PVdim,Kedim);
        E=zeros(PVdim,PVdim);
        for j=1:npts_c % loop over all deviatoric-strain quadrature points
            c =Ns_c{j}'*x;% physical location of the quadrature point
            J = x' * Nders_c{j};% We compute the Jacobian matrix
            if (~Rm_constant)% do I need to evaluate the local material orientation?
                if (~isempty(labels )),  Rm =Rmh(c,J,labels(i));
                else,                    Rm =Rmh(c,J,[]);                end
            end
            Ndersp = Nders_c{j}/J;% derivatives wrt global coor
            Jac = Jacobian_volume(fes,conn, Ns_c{j}, J, x);
            if (Jac<=0),error('Non-positive Jacobian');end
            if (~D_constant)
                D = tangent_moduli(mat,struct('xyz',c));
                [m1,Id,Iv]=constrained_mat_data(self,D);
            end
            B = self.hBlmat(self,Ns_c{j},Ndersp*Rm,c,Rm);% strain-displacement
            CT = CT + (Ns_pv_c{j}*(Jac*w_c(j)))*m1'*B;
            E = E + Ns_pv_c{j}*(Ns_pv_c{j}*(Jac*w_c(j)))';
        end
        W=E\CT;
        for j=1:npts_u % loop over all volumetric-strain quadrature points
            c =Ns_u{j}'*x;% physical location of the quadrature point
            u_c = transpose(Ns_u{j})*Us(conn,:);
            J = x' * Nders_u{j};% We compute the Jacobian matrix
            if (~Rm_constant)% do I need to evaluate the local material orientation?
                if (~isempty(labels )),  Rm =Rmh(c,J,labels(i));
                else,                    Rm =Rmh(c,J,[]);                end
            end
            if (~outputRm_constant)% do I need to evaluate the output orientation?
                if (~isempty(labels )),  outputRm =outputRmh(c,J,labels(i));
                else,                    outputRm =outputRmh(c,J,[]);    end
            end
            Ndersp = Nders_u{j}/J;% derivatives wrt global coor
            Jac = Jacobian_volume(fes,conn, Ns_u{j}, J, x);
            if (Jac<=0),error('Non-positive Jacobian');end
            if (~D_constant)
                D = tangent_moduli(mat,struct('xyz',c));
                [m1,Id,Iv]=constrained_mat_data(self,D);
            end
            B = self.hBlmat(self,Ns_u{j},Ndersp*Rm,c,Rm);% strain-displacement
            Bbar=Id*B +  (1/3)*m1*Ns_pv_u{j}'*W;
            context.strain = Bbar*U;% strain in material coordinates
            context.dT = transpose(Ns_u{j})*dT;
            context.xyz =c;
            [out,ignore] = update(mat, [], context);
            switch context.output
                case 'Cauchy'
                    out =mat.stress_vector_rotation(Rm')*out;%  To global coordinate system
                    if (~outputRm_identity)
                        if (~outputRm_constant)
                            c =mean(X);% physical location of the quadrature point
                            if (~isempty(labels )),  outputRm=outputRmh(c,[],labels(i));%  No Jacobian matrix?
                            else,                    outputRm =outputRmh(c,[],[]);                end
                        end
                        out =mat.stress_vector_rotation(outputRm)*out;% To output coordinate system
                    end
            end
            if ~isempty (inspector)
                idat =feval(inspector,idat,out,c,u_c,pc_u(j,:));
            end 
        end
    end
    return;
end
