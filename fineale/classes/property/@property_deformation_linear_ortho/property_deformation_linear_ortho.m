classdef property_deformation_linear_ortho < property_base
    % Class of properties for linearly elastic orthotropic materials.
    %
    %
    
    properties
        % Young's moduli
        E1=[];      E2=[];      E3=[];
        % Poisson ratio
        nu12=[];    nu13=[];    nu23=[];
        % Shear moduli
        G12=[];     G13=[];     G23=[];
        % Thermal expansion coefficients
        alpha1= []; alpha2= []; alpha3= [];
    end
    
    properties (Hidden, SetAccess = private)
        Precomputed_D = [];
    end
    
    methods
        
        function self = property_deformation_linear_ortho (Parameters)
            % Constructor.
            % Parameters:
            % required fields
            %     E1, E2, E3=Young's modulus
            %     nu12, nu13, nu23=Poisson ratio
            %     G12,  G13, G23=Shear modulus
            % optional fields
            %     alpha = coefficient of thermal expansion, default 0.0
            %     nconstrained= number of constrained modes (default is 1)
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            
            if (nargin < 1)
                Parameters =struct( [] );
            end
            self = self@property_base(Parameters);
            if ( ~isfield( Parameters, 'E1'))
                return
            end
            E1=Parameters.E1;
            E2=Parameters.E2;
            E3=Parameters.E3;
            G12=Parameters.G12;
            G13=Parameters.G13;
            G23=Parameters.G23;
            nu12=Parameters.nu12;
            nu13=Parameters.nu13;
            nu23=Parameters.nu23;
            % compute the compliance, then the stiffness
            compliance =[1/E1      -nu12/E1    -nu13/E1  0   0   0;...
                -nu12/E1     1/E2      -nu23/E2  0   0   0;...
                -nu13/E1   -nu23/E2       1/E3   0   0   0;...
                0           0           0 1/G12 0   0;...
                0           0           0   0 1/G13 0;...
                0           0           0   0   0 1/G23];
            if rank(compliance)<6
                error(' Singular compliance?')
            end
            if prod(eig (compliance)) <= 0
                error(' Non positive-definite compliance?')
            end
            self.Precomputed_D= inv(compliance);
            self.E1=E1;self.E2=E2;self.E3=E3;
            self.nu12=nu12;self.nu13=nu13;self.nu23=nu23;
            self.G12=G12;self.G13=G13;self.G23=G23;
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
