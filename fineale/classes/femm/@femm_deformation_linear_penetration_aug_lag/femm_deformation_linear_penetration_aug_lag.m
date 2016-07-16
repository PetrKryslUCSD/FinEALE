% Class for  penalty and Lagrangean-multiplier forces associated with penetration of contact
% surface
% 
%
classdef femm_deformation_linear_penetration_aug_lag < femm_deformation_linear
    
    properties %(Hidden, SetAccess = private)
        surface_data = [];% Whatever you wish to have passed into the get_penetration() function.
        get_penetration = [];% Function handle: function signature
        %  [penetration,normal] = get_penetration(surface_data,qpx,qpu);
        % where qpx,qpu are the location and the displacement at the quadrature point.
        penalty =  1.0;% Penalty parameter for the contact force: Force per unit area and unit penetration.
        lm=[];% Lagrange-multiplier forces, one per quadrature point
    end
    
    methods % constructor
        
        function self = femm_deformation_linear_penetration_aug_lag (Parameters)
            if nargin <1
                Parameters = struct('surface_data',[]);
            end
            self = self@femm_deformation_linear(Parameters);
            if (isfield(Parameters,'get_penetration'))
                self.get_penetration = Parameters.get_penetration;
                self.penalty =  Parameters.penalty;
            end
            self.surface_data = [];
            if (isfield(Parameters,'surface_data'))
                self.surface_data = Parameters.surface_data;
            end
            [npts Ns Nders w] = integration_data (self);
            self.lm= zeros(count( self.fes),npts);       
        end
        
        function [F, updatedself] = contact_loads(self, assembler, geom, u)
            % Compute the system load vector due to penetration of contact surface.
            %
            % function F = contact_loads(self, assembler, geom, u)
            %
            %     geom=geometry field
            %     u=displacement field
            %
            fes = self.fes;% grab the finite elements to work on
            % Integration rule: compute the data needed for  numerical quadrature
            [npts Ns Nders w] = integration_data (self);
            conns = fes.conn; % connectivity
            labels = fes.label; % connectivity
            xs =geom.values;
            Us=u.values;
            updatedself=self;
            lm= self.lm;
            %             Precompute for efficiency
            dim=u.dim; nfens=fes.nfens; Id =eye(dim); NexpT={};
            for    qp =1:npts
                Nexp=zeros(dim,nfens);
                for l = 1:nfens
                    Nexp(1:dim,(l-1)*dim+1:(l)*dim)=Id*Ns{qp}(l);
                end;
                NexpT{qp}=Nexp';
            end; clear qp
            % Prepare assembler
            Kedim =u.dim*fes.nfens;
            start_assembly(assembler, u.nfreedofs);
            % Now loop over all Finite elements
            for i=1:size(conns,1)
                conn =conns(i,:);
                dofnums =reshape(u,gather_dofnums(u,conn));
                X=xs(conn,:);
                U=Us(conn,:);
                Fe = zeros(Kedim,1);       % element force vector
                for    qp =1:npts
                    J = Jacobian_matrix(fes,Nders{qp},X);
                    Jac = Jacobian_surface(fes,conn, Ns{qp}, J, X);
                    qpx=Ns{qp}'*X;
                    qpu=Ns{qp}'*U;
                    [penetration,normal] = self.get_penetration(self.surface_data,qpx,qpu);
                    % Predict the Lagrange multiplier
                    lmpred=lm(i,qp)+(penetration)*Jac*w(qp)*self.penalty;
                    if (lmpred<0) % Macaulay-bracket the Lagrange multiplier
                        lmpred=0.0;
                    end
                    if (lmpred>0)% if there is any contribution to the force, add it now
                        Fe =  Fe + NexpT{qp}*(normal'*lmpred);
                    end
                    lm(i,qp)=lmpred;
                end; clear qp
                assemble(assembler, Fe, dofnums);
            end
            F = make_vector (assembler);
            updatedself.lm=lm;
        end
        
        function K = stiffness (self, assembler, geom, u)
            % Compute the stiffness matrix corresponding to the contact penalty.
            %
            % function K = stiffness (self, assembler, geom, u)
            %
            %     geom=reference geometry field
            %     u=displacement field
            %    Returns the system stiffness matrix.
            %
            fes = self.fes;% grab the finite elements to work on
            % Integration rule: compute the data needed for  numerical quadrature
            [npts Ns Nders w] = integration_data (self);
            conns = fes.conn; % connectivity
            labels = fes.label; % connectivity
            xs =geom.values;
            Us=u.values;
            lm= self.lm;
            %             Precompute for efficiency
            dim=u.dim; nfens=fes.nfens; Id =eye(dim); NexpS={}; NexpT={};
            for    qp =1:npts
                Nexp=zeros(dim,nfens);
                for l = 1:nfens
                    Nexp(1:dim,(l-1)*dim+1:(l)*dim)=Id*Ns{qp}(l);
                end;
                NexpS{qp}=Nexp;
                NexpT{qp}=Nexp';
            end; clear qp
            % Prepare assembler
            Kedim =u.dim*fes.nfens;
            start_assembly(assembler, Kedim, Kedim, size(conns,1), u.nfreedofs, u.nfreedofs);
            % Now loop over all fes in the block
            for i=1:size(conns,1)
                conn =conns(i,:);
                dofnums =reshape(u,gather_dofnums(u,conn));
                X=xs(conn,:);
                U=Us(conn,:);
                Ke = zeros(Kedim);       % element stiffness matrix
                for    qp =1:npts
                    J = Jacobian_matrix(fes,Nders{qp},X);
                    Jac = Jacobian_surface(fes,conn, Ns{qp}, J, X);
                    qpx=Ns{qp}'*X;
                    qpu=Ns{qp}'*U;
                    [penetration,normal] = self.get_penetration(self.surface_data,qpx,qpu);
                    % Predict the Lagrange multiplier
                    lmcur=lm(i,qp)+(penetration)*Jac*w(qp)*self.penalty;
                    if (lmcur>0)
                        Ke = Ke + NexpT{qp}*(self.penalty*Jac*w(qp)*normal'*normal)*NexpS{qp};
                    end
                end; clear qp
                assemble_symmetric(assembler, Ke, dofnums);
            end
            K = make_matrix (assembler);
        end
        
    end
    
end