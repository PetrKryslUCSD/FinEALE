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
    [npts_v Ns_v Nders_v w_v pc_v] = integration_data_constrained(self);
    [npts_s Ns_s Nders_s w_s pc_s] = integration_data_unconstrained(self);
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
    ms=[];
    % Retrieve data for efficiency
    conns = fes.conn; % connectivity
    labels = fes.label; % finite element labels
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
        for j=1:npts_v % loop over all volumetric-strain quadrature points
            c =Ns_v{j}'*x;% physical location of the quadrature point
            u_c = transpose(Ns_v{j})*Us(conn,:);
            J = x' * Nders_v{j};% We compute the Jacobian matrix
            if (~Rm_constant)% do I need to evaluate the local material orientation?
                if (~isempty(labels )),  Rm =Rmh(c,J,labels(i));
                else,                    Rm =Rmh(c,J,[]);                end
            end
            if (~outputRm_constant)% do I need to evaluate the output orientation?
                if (~isempty(labels )),  outputRm =outputRmh(c,J,labels(i));
                else,                    outputRm =outputRmh(c,J,[]);    end
            end
            Ndersp = Nders_v{j}/J;% derivatives wrt global coor
            Jac = Jacobian_volume(fes,conn, Ns_v{j}, J, x);
            if (Jac<=0),error('Non-positive Jacobian');end
            B=self.hBlmat(self,Ns_v{j},Ndersp*Rm,c,Rm);%  strains in material coordinates, displacements in global coordinates
            context.strain = B*U;%strain in the  material coordinates
            context.dT = transpose(Ns_v{j})*dT;
            context.xyz =c;
            context.Rm =Rm;
            [out,ignore] = update(mat, [], context);
            switch context.output
                case 'Cauchy'
                    out =mat.stress_vector_rotation(Rm')*out;%  To global coordinate system
                    if (~outputRm_identity)
                        out =mat.stress_vector_rotation(outputRm)*out;% To output coordinate system
                    end
            end
            if ~isempty (inspector)
                idat =feval(inspector,idat,out,c,u_c,pc_v(j,:));
            end 
        end
    end
    return;
end

