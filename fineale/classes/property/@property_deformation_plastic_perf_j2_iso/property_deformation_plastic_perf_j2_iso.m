% Class of properties for linearly elastic isotropic materials. 
%
classdef property_deformation_plastic_perf_j2_iso < property_deformation_linear_iso
    
    properties
        sigma_y =  []; % Yield stress
    end
    
    methods
        
        function self = property_deformation_plastic_perf_j2_iso (Parameters)
            % Constructor.
            % Parameters:
            % Required fields
            %     E=Young's modulus
            %     nu=Poisson ratio
            % Optionally one may supply the shear modulus
            %     G=shear modulus
            %       The Poisson ratio will be then computed from Young's modulus and
            %       the shear modulus
            % Optional fields
            %     alpha = coefficient of thermal expansion, default 0.0
            %
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if (nargin < 1)
                Parameters =struct( [] );
            end
            self = self@property_deformation_linear_iso(Parameters);
            if ( ~isfield( Parameters, 'sigma_y'))
                return
            end
            self.sigma_y = 0.0;
            if isfield(Parameters,'sigma_y')
                self.sigma_y = Parameters.sigma_y;
                if (self.sigma_y <= 0)
                    error('Non-positive sigma_y!');
                end
            end
            
        end
        
    end

end