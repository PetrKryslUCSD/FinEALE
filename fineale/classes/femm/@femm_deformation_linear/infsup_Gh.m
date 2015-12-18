% Compute the pressure norm for the inf-sup numerical test.
%
% function Gh = infsup_Gh (self, assembler, geom, u)
%
function Gh = infsup_Gh (self, assembler, geom, u)
fes = self.fes;% grab the finite elements to work on
    % Integration rule: compute the data needed for  numerical quadrature
    [npts Ns Nders w] = integration_data (self);
    % Material orientation
    Rm_constant = is_material_orientation_constant(self);% if not constant, need to compute  at each point
    if (~Rm_constant)
        Rmh = self.Rm;% handle to a function  to evaluate Rm
    else
        Rm = self.Rm;% constant material orientation matrix
    end    
    % Material
    mat = self.material;
    D_constant = are_tangent_moduli_constant (mat);
    if (D_constant)
        D = tangent_moduli(mat,[]);%  Moduli in material orientation
    end
    % Retrieve data for efficiency
    conns = fes.conn; % connectivity
    labels = fes.label; % connectivity
    xs =geom.values;
    % Prepare assembler
    Kedim =u.dim*fes.nfens;
    start_assembly(assembler, Kedim, Kedim, size(conns,1), u.nfreedofs, u.nfreedofs);
    % Now loop over all fes in the block
    for i=1:size(conns,1)
        conn =conns(i,:);
        x=xs(conn,:);
        dofnums =reshape(u,gather_dofnums(u,conn));
        Ke =zeros(Kedim);;
        for j=1:npts
            c =Ns{j}'*x;% physical location of the quadrature point
            J = x' * Nders{j};% We compute the Jacobian matrix
            if (~Rm_constant)% do I need to evaluate the local material orientation?
                if (~isempty(labels )),  Rm =Rmh(c,J,labels(i));
                else,                    Rm =Rmh(c,J,[]);                end
            end
            Ndersp = Nders{j}/J;% derivatives wrt global coor
            Jac = Jacobian_volume(fes,conn, Ns{j}, J, x);
            if (Jac<=0),error('Non-positive Jacobian');end
            B = self.hdivmat(self,Ns{j},Ndersp*Rm,c,Rm);% strain-displacement
            Ke = Ke + B'*B * (Jac * w(j));
        end
        assemble_symmetric(assembler, Ke, dofnums);
    end
    Gh = make_matrix (assembler);
end

% Due to  strain-displacement matrices below are equivalent
%             B = self.hBlmat(self,Ns{j},Ndersp,c,[]);% strain-displacement
%             strain_vector_rotation(mat,Rm')*B-self.hBlmat(self,Ns{j},Ndersp*Rm,c,Rm)

%
%     gcells = get(self.feblock,'gcells');
%     nfens = get(gcells,'nfens');
%     % Integration rule
%     integration_rule = get(self.feblock, 'integration_rule');
%     pc = get(integration_rule, 'param_coords');
%     w  = get(integration_rule, 'weights');
%     npts_per_gcell = get(integration_rule, 'npts');
%     for j=1:npts_per_gcell
%         Ns{j} = bfun(gcells,pc(j,:));
%         Nders{j} = bfundpar(gcells,pc(j,:));
%     end
%     % Material orientation /directions
%     Rm = get(self,'Rm');
%     eval_Rm = 0;
%     if strcmp(class(Rm),'function_handle')
%         eval_Rm = true; Rmh =Rm;
%     end
%     % Material
%     mat = get(self.feblock, 'mater');
%     % Now loop over all gcells in the block
%     conns = get(gcells, 'conn'); % connectivity
%     Ke = cell(size(conns,1),1);
%     eqnums = cell(size(conns,1),1);
%     for i=1:size(conns,1)
%         eqnums{i} =gather(u,conns(i,:),'eqnums');
%         Ke{i} =zeros(get(geom,'dim')*nfens);
%     end
%     xs =gather(geom,(1:get(geom,'nfens')),'values','noreshape');
%     for i=1:count(gcells)
%         conn =conns(i,:);
%         x=xs(conn,:);
%         % Loop over all integration points
%         for j=1:npts_per_gcell
%             c =Ns{j}'*x;% physical location of the quadrature point
%             J = x' * Nders{j};% We compute the Jacobian matrix
%             if (eval_Rm)% do I need to evaluate the local material orientation?
%                 Rm =Rmh(c,J);
%             end
%             if (isempty(Rm))
%                 Ndersp = Nders{j}/J;% derivatives wrt local coor
%             else
%                 Ndersp = Nders{j}/(Rm'*J);% derivatives wrt local coor
%             end
%             Jac = Jacobian_volume(gcells,conn, Ns{j}, J, x);
%             if (Jac<=0),error('Non-positive Jacobian');end
%             B = self.hdivmat(self,Ns{j},Ndersp,c,Rm);% strain-displacement
%             Ke{i} = Ke{i} + B'*B * (Jac * w(j));
%         end
%     end
%     ems  = elematset(struct('mat',{Ke}, 'eqnums',{eqnums}));
% end