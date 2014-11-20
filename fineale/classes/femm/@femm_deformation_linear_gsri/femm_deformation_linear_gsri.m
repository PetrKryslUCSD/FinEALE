classdef femm_deformation_linear_gsri < femm_deformation_linear
% Class for the small displacement, small strain deformation model with 
% generalized selective reduced integration (GSRI).
%
% Reference: Generalized Selective Reduced Integration and B-bar Finite Element Methods for Anisotropic Elasticity,
% P. Krysl, J. Nov\'ak, S. Oberrecht, International Journal for numerical methods in engineering, submitted 2013.
    
    properties
        fraction_constrained=1;% How much of the constrained energy should be included?
        integration_rule_constrained= [];% Integration rule for constrained terms
        fraction_unconstrained=1;% How much of the unconstrained energy should be included?
        integration_rule_unconstrained= [];% Integration rule for unconstrained terms
        nconstrained = 1;% how many tangent-moduli modes should be considered constrained?
    end
     
    methods % constructor
        
        function self = femm_deformation_linear_gsri (Parameters)
            % Constructor.
            % Parameters:
            % Options: those recognized by model_deformation_linear plus
            %    N = Outer normal, either as an array of numbers,or as a function handle.
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin <1
                Parameters = [];
            end
            self = self@femm_deformation_linear(Parameters);
            if isfield(Parameters,'fraction_constrained')
                self.fraction_constrained = Parameters.fraction_constrained;
            end
            self.integration_rule_constrained = self.integration_rule;
            if isfield(Parameters,'integration_rule_constrained')
                self.integration_rule_constrained = Parameters.integration_rule_constrained;
            end
            if isfield(Parameters,'fraction_unconstrained')
                self.fraction_unconstrained = Parameters.fraction_unconstrained;
            end
            self.integration_rule_unconstrained = self.integration_rule;
            if isfield(Parameters,'integration_rule_unconstrained')
                self.integration_rule_unconstrained = Parameters.integration_rule_unconstrained;
            end
            self.nconstrained = 1;
            if isfield(Parameters,'nconstrained')
                self.nconstrained = Parameters.nconstrained;
            end
        end
        
    end
    
    methods (Hidden)
         
         
        function [D_constrained,D_unconstrained]=constrained_mat_data(self,D)
            % Calculate the constrained directions.
            %
            % function [m1]=constrained_directions(self,D)
            %
            % Arguments
            %    self=property
            [DV,DL]=eig(D);
            [ignore, iix] =sort(diag(DL));
            m1 = sqrt(3)*DV(:,iix(end-self.nconstrained+1:end));
            Iv =1/3*(m1 * m1');
            Id =eye(6)-Iv;
            D_constrained=Iv*D*Iv;
            D_unconstrained=Id*D*Id;
        end
        
        function [npts Ns gradNparams w pc] = integration_data_constrained (self)
        % Calculate the data needed for  numerical quadrature. Constrained terms.
            pc = self.integration_rule_constrained.param_coords;
            w  =  self.integration_rule_constrained.weights ;
            npts = self.integration_rule_constrained.npts;
            % Precompute basis f. values + basis f. gradients wrt parametric coor
            clear Ns Nders
            for j=1:npts
                Ns{j} = bfun(self.fes,pc(j,:));
                gradNparams{j} = bfundpar(self.fes,pc(j,:));
            end
        end
        
        function [npts Ns gradNparams w pc] = integration_data_unconstrained (self)
        % Calculate the data needed for  numerical quadrature. Unconstrained terms.
            pc = self.integration_rule_unconstrained.param_coords;
            w  =  self.integration_rule_unconstrained.weights ;
            npts = self.integration_rule_unconstrained.npts;
            % Precompute basis f. values + basis f. gradients wrt parametric coor
            clear Ns Nders
            for j=1:npts
                Ns{j} = bfun(self.fes,pc(j,:));
                gradNparams{j} = bfundpar(self.fes,pc(j,:));
            end
        end
        
    end
    
end
