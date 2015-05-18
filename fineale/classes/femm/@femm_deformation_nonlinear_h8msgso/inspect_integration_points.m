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
    end
end
if ~isfield(context,'output')
    context.output=  'Cauchy';% default output
end
% Retrieve data for efficiency
conns = fes.conn; % connectivity
labels = fes.label; % finite element labels
Xs =geom.values;
Us =u.values;
context.F= [];
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
    U=Us(conn,:);
    x1 = X + U; % current coordinates
    dT =dTs(conn,:);
    % First we calculate  the mean basis function gradient matrix and the volume of the element
    gradN_mean =zeros(self.fes.nfens,u.dim); V=0;
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
    F1bar =x1'*gradN_mean;% wrt global material coordinates
    context.F =Rm'*F1bar*Rm;% Deformation gradient wrt  material orientation 
    context.dT = transpose(Ns{j})*dT;
    context.xyz =mean(X);
    [out,ignore] = update(self.material, self.matstates{i}, context);
    switch context.output
        case 'Cauchy'
            out =self.material.stress_vector_rotation(Rm')*out;%  To global coordinate system
            if (~outputRm_identity)
                out =self.material.stress_vector_rotation(outputRm)*out;% To output coordinate system
            end
    end
    if ~isempty (inspector)
        idat =feval(inspector,idat,out,mean(X),transpose(Ns{j})*U,[0,0,0]);
    end
end
return;
end





