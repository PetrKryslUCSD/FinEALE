% Class of properties for hyper elastic neo-Hookean material. 
%
classdef property_deformation_neohookean < property_base
    
    properties
        E = [];% Young's modulus
        nu = [];% Poisson ratio
        alpha = [];% Coefficient of thermal expansion
    end
    
    properties (Hidden, SetAccess = private)
            lambda = [];
            mu     = [];
    end
    
    methods
        
        function self = property_deformation_neohookean (Parameters)
            % Constructor.
            % Parameters:
            % Required fields
            %     E=Young's modulus
            %     nu=Poisson ratio
            % Optional fields
            %     alpha = coefficient of thermal expansion, default 0.0
            %
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if (nargin < 1)
                Parameters =struct( [] );
            end
            self = self@property_base(Parameters);
            if ( ~isfield( Parameters, 'E'))
                return
            end
            self.E = Parameters.E;
            self.nu =  0;
            if isfield(Parameters,'nu')
                self.nu = Parameters.nu;
                if (isnumeric(self.nu))
                    if (self.nu < 0)
                        error('Negative Poisson ratio!');
                    elseif (self.nu >= 0.5)
                        error('Incompressible material!');
                    end
                end
            end
            self.alpha = 0.0;
            if isfield(Parameters,'alpha')
                self.alpha= Parameters.alpha;
            end
            % precompute
            %             self.lambda = E * nu / (1 + nu) / (1 - 2*(nu));
            %             self.mu     = E / (2 * (1 + nu));
        end
        
        function val = thermal_expansion_coefficients (self)
        % Return the coefficients of thermal expansion in all three directions.
            val = self.alpha*ones(3, 1);
        end
        
        function val = are_tangent_moduli_constant (self)
        % Is the material stiffness matrix independent of location?
        %
        % Is the material stiffness matrix independent of location (constant,
        % corresponding to a homogeneous material)?
            val = true;
        end
        
        
    end
    
    
end
