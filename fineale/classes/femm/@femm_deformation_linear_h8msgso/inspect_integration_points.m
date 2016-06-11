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
[npts, Ns, gradNparams, w, pc] = integration_data (self);;
gradN =cell(8,1); Jac=cell(8,1);
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
    else
        outputRm = Rm;
    end
end
if ~isfield(context,'output')
    context.output=  'Cauchy';% default output
end
% Retrieve data for efficiency
conns = fes.conn; % connectivity
labels = fes.label; % finite element labels
Xs =geom.values;
Un1s =un1.values;
context.Fn1= [];
if isempty(dT)
    dTs=zeros(geom.nfens,1);
else
    dTs=dT.values;
end
% Now loop over selected fes
for m=1:length(fe_list)
    i=fe_list(m);
    conn = conns(i,:); % connectivity
    X=Xs(conn,:);
    U=reshape(un1, gather_values(un1,conn));
    dT =dTs(conn,:);
    % First we calculate  the mean basis function gradient matrix and the volume of the element
    gradN_mean =zeros(self.fes.nfens,un1.dim); V=0;
    for j=1:npts % loop over all quadrature points
        J = X' * gradNparams{j};% Jacobian matrix wrt reference coordinates
        gradN{j} = gradNparams{j}/J;% derivatives wrt reference coor
        Jac{j} = det(J);        if (Jac{j}<=0),error('Non-positive Jacobian');end
        dV=(Jac{j}*w(j));
        gradN_mean =gradN_mean+gradN{j}*dV;
        V =V+dV;
    end
    gradN_mean=gradN_mean/V; % Finish the calculation of the mean basis function gradient matrix
    %     Material orientation?
    if (~Rm_constant)% do I need to evaluate the local material orientation?
        c =mean(X);% physical location of the quadrature point
        if (~isempty(labels )),  Rm =Rmh(c,[],labels(i));%  No Jacobian matrix?
        else,                    Rm =Rmh(c,[],[]);                end
    end
    % Now we calculate the mean deformation gradient dx/dX   
    Bbar = self.hBlmat(self,[],gradN_mean,[],[]);% strain-displacement d/dX
    context.strain =Bbar*U;% Strain wrt  material orientation 
    context.dT = transpose(Ns{j})*dT;
    context.xyz =mean(X);
    out = state(self.material, [], context);
    switch context.output
        case 'Cauchy'
            out =self.material.stress_vector_rotation(Rm')*out;%  To global coordinate system
            if (~outputRm_identity)
                if (~outputRm_constant)
                    c =mean(X);% physical location of the quadrature point
                    if (~isempty(labels )),  outputRm=outputRmh(c,[],labels(i));%  No Jacobian matrix?
                    else,                    outputRm =outputRmh(c,[],[]);                end
                end
                out =self.material.stress_vector_rotation(outputRm)*out;% To output coordinate system
            end
    end
    if ~isempty (inspector)
        idat =feval(inspector,idat,out,mean(X),transpose(Ns{j})*Un1s(conn,:),[0,0,0]);
    end
end
return;
end





