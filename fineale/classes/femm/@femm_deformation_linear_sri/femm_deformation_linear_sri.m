classdef femm_deformation_linear_sri < femm_deformation_linear
% Class for the small displacement, small strain deformation model
% with selective reduced integration.
%
%
    
    properties
        integration_rule_volumetric= [];% integration rule for volumetric terms
        integration_rule_deviatoric= [];% integration rule for deviatoric terms
        split ='lambda';% kind of split of energy
    end
    
    
    methods % constructor
        
        function self = femm_deformation_linear_sri (Parameters)
            % Constructor.
            % Parameters:
            % Options: those recognized by model_deformation_linear plus
            %    N = Outer normal, either as an array of numbers,or as a function handle.
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin <1
                Parameters = [];
            end
            self = self@femm_deformation_linear(Parameters);
            self.integration_rule_volumetric = self.integration_rule;
            if isfield(Parameters,'integration_rule_volumetric')
                self.integration_rule_volumetric = Parameters.integration_rule_volumetric;
            end
            self.integration_rule_deviatoric = self.integration_rule;
            if isfield(Parameters,'integration_rule_deviatoric')
                self.integration_rule_deviatoric = Parameters.integration_rule_deviatoric;
            end
            if isfield(Parameters,'split')
                self.split = Parameters.split;
            end
            if ~(strcmp(self.split,'lambda') ||  strcmp(self.split,'bulk'))
                error ('Unknown split')
            end
            
        end
        
    end
    
    methods (Hidden)
        
        
        function [npts Ns gradNparams w pc] = integration_data_volumetric (self)
            % Calculate the data needed for  numerical quadrature.
            %
            % function [npts Ns gradNparams w pc] = integration_data_volumetric (self)
            pc = self.integration_rule_volumetric.param_coords;
            w  =  self.integration_rule_volumetric.weights ;
            npts = self.integration_rule_volumetric.npts;
            % Precompute basis f. values + basis f. gradients wrt parametric coor
            clear Ns Nders
            for j=1:npts
                Ns{j} = bfun(self.fes,pc(j,:));
                gradNparams{j} = bfundpar(self.fes,pc(j,:));
            end
        end
        
        function [npts Ns gradNparams w pc] = integration_data_deviatoric (self)
            % Calculate the data needed for  numerical quadrature.
            %
            % function [npts Ns gradNparams w pc] = integration_data_deviatoric (self)
            %
            pc = self.integration_rule_deviatoric.param_coords;
            w  =  self.integration_rule_deviatoric.weights ;
            npts = self.integration_rule_deviatoric.npts;
            % Precompute basis f. values + basis f. gradients wrt parametric coor
            clear Ns Nders
            for j=1:npts
                Ns{j} = bfun(self.fes,pc(j,:));
                gradNparams{j} = bfundpar(self.fes,pc(j,:));
            end
        end
        
        
    end
    
end
