classdef property_deformation_linear_aniso < property_base
    % Class of properties for linearly elastic general anisotropic materials.
    %
    %
    
    properties
    % Material stiffness matrix
        D=[];
        % Thermal expansion coefficients
        alpha1= []; alpha2= []; alpha3= [];
    end
    
    properties (Hidden, SetAccess = private)
        Precomputed_D = [];
    end
    
    methods
        
        function self = property_deformation_linear_aniso (Parameters)
            % Constructor.
            % Parameters:
            % required fields
            %     D=Material stiffness matrix
            % optional fields
            %     alpha = coefficient of thermal expansion, default 0.0
            %     nconstrained= number of constrained modes (default is 1)
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            
            if (nargin < 1)
                Parameters =struct( [] );
            end
            self = self@property_base(Parameters);
            if ( ~isfield( Parameters, 'D'))
                return
            end
            D=Parameters.D;
            self.Precomputed_D= D;
            % thermal properties
            self.alpha1=0;
            self.alpha2=0;
            self.alpha3=0;
            if isfield(Parameters,'alpha1')
                self.alpha1= Parameters.alpha1;
                self.alpha2= Parameters.alpha2;
                self.alpha3= Parameters.alpha3;
            end
        end
        
        function val = thermal_expansion_coefficients (self)
            % Return the coefficients of thermal expansion in all three directions.
            val = [self.alpha1;self.alpha2;self.alpha3];
        end
        
        function val = are_tangent_moduli_constant (self)
            % Is the material stiffness matrix independent of location?
            % Is the material stiffness matrix independent of location (constant,
            % corresponding to a homogeneous material)?
            val = true;
        end
        
        
        function val = tangent_moduli(self, context)
            % Calculate the material stiffness matrix.
            %
            % function val = tangent_moduli(self, context)
            %
            % Arguments
            %    self=property
            %    context=structure with optional field:
            %                kind = string with possible values
            %                       constrained, unconstrained
            %
            % Output
            %     D=matrix 6x6 in the local material orientation Cartesian basis
            %
            if isfield(context,'kind')
                switch context.kind
                    otherwise
                        val = [];
                end
            else
                val= self.Precomputed_D;
            end
            return;
        end
        
        
    end
    
end
