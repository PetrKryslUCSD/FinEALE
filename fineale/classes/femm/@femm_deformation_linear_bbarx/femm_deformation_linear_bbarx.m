classdef femm_deformation_linear_bbarx < femm_deformation_linear
% Class for the small displacement, small strain deformation model with 
% generalized B-bar formulation with stabilization. 
% Reference:  will
% B-bar FEMs for anisotropic elasticity
% S.P. Oberrecht, J. Novák andP. Krysl,
% International Journal for Numerical Methods in Engineering
% Volume 98, Issue 2, pages 92–104, 13 April 2014 

    properties
        integration_rule_constrained= [];% Integration rule for constrained terms (the E,C matrices).
        integration_rule_unconstrained= [];% Integration rule for unconstrained terms (the  stiffness matrix itself).
        % Function handle, function to calculate the values of pressure/volumetric strain basis functions at quadrature points.
          % The function could be supplied for instance like this
          % @(p)[1;p(1);p(2);p(3);]
           % which would be complete linear function in parametric coordinates.
        pv_bfun=[];
        nconstrained = 1;% how many tangent-moduli modes should be considered constrained?
        psi=[]; %  Parameter  in the split of the material stiffness moduli.
    end
     
    properties (Access = private)
        Psi_mat=[]; %  Parameter  in the split of the material stiffness moduli.
    end
 
    methods % constructor
        
        function self = femm_deformation_linear_bbarx (Parameters)
            % Constructor.
            % Parameters:
            % Options: those recognized by model_deformation_linear plus
            %    N = Outer normal, either as an array of numbers,or as a function handle.
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin <1
                Parameters = [];
            end
            self = self@femm_deformation_linear(Parameters);
            self.integration_rule_constrained = self.integration_rule;% Default: take it from the parent
            if isfield(Parameters,'integration_rule_constrained')
                self.integration_rule_constrained = Parameters.integration_rule_constrained;
            end
            self.integration_rule_unconstrained = self.integration_rule;% Default: take it from the parent
            if isfield(Parameters,'integration_rule_unconstrained')
                self.integration_rule_unconstrained = Parameters.integration_rule_unconstrained;
            end
            self.pv_bfun = @(Parametric_coordinate)1;% Default pressure/volumetric strain basis functions.
            if isfield(Parameters,'pv_bfun')
                self.pv_bfun = Parameters.pv_bfun;
            end
            self.nconstrained = 1;
            if isfield(Parameters,'nconstrained')
                self.nconstrained = Parameters.nconstrained;
            end
            if (isempty(self.pv_bfun))
                self.pv_bfun=@(p)polynomial_basis(3,0,p)*eye(self.nconstrained);
            end
            self.psi = [];
            if isfield(Parameters,'psi')
                self.psi = Parameters.psi;
            end
        end
        
    end
    
    methods (Hidden)
         
         
        function [m1,Id,Psi_mat]=constrained_mat_data(self,D)
            % Calculate the constrained directions.
            %
            % function [m1]=constrained_directions(self,D)
            %
            % Arguments
            %    self=property
            if (~isempty(self.nconstrained))
                [DV,DL]=eig(D);
                [DLa, iix] =sort(diag(real(DL)));
                cDLa=DLa(end-self.nconstrained+1:end);
                fDLa=DLa(1:end-self.nconstrained);
                m1 = sqrt(3)*DV(:,iix(end-self.nconstrained+1:end));
                if (isempty(self.psi))
                    psi=1-min(fDLa)/max(cDLa);
                else
                    psi= self.psi;
                end
                Psi_mat=diag(1-psi*(min(DLa)./cDLa));
                Psi_mat(end,end)=1.0; % Totally wipe out the stiffest mode
            else
                m1=[1,1,1,0,0,0]';
                Psi_mat=eye(size(m1,2));
            end
            Id =eye(6)- m1 * 1/3* (eye(size(m1,2))-diag(sqrt(1-diag(Psi_mat)))) * m1';
            %                         fD=eig(Id*D*Id);
            %                         semilogy(fD,'b-o');
            %                         hold on
            %                         cD=eig(D-Id*D*Id);
            %                         semilogy(cD,'r-s');
        end
        
         function [npts Ns gradNparams Ns_pv w pc] = integration_data_constrained (self)
        % Calculate the data needed for  numerical quadrature.
        % The matrices E,C are calculated using this quadrature.
           pc = self.integration_rule_constrained.param_coords;
            w  =  self.integration_rule_constrained.weights ;
            npts = self.integration_rule_constrained.npts;
            % Precompute basis f. values + basis f. gradients wrt parametric coor
            clear Ns Nders
            for j=1:npts
                Ns{j} = bfun(self.fes,pc(j,:));
                gradNparams{j} = bfundpar(self.fes,pc(j,:));
                Ns_pv{j} = self.pv_bfun(pc(j,:));
            end
        end
        
        function [npts Ns gradNparams Ns_pv w pc] = integration_data_unconstrained (self)
        % Calculate the data needed for  numerical quadrature.
        % The stiffness matrix is calculated using this quadrature.
            pc = self.integration_rule_unconstrained.param_coords;
            w  =  self.integration_rule_unconstrained.weights ;
            npts = self.integration_rule_unconstrained.npts;
            % Precompute basis f. values + basis f. gradients wrt parametric coor
            clear Ns Nders
            for j=1:npts
                Ns{j} = bfun(self.fes,pc(j,:));
                gradNparams{j} = bfundpar(self.fes,pc(j,:));
                Ns_pv{j} = self.pv_bfun(pc(j,:));
            end
        end


    end
    
end
