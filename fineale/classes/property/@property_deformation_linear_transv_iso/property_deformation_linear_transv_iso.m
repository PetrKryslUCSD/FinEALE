classdef property_deformation_linear_transv_iso < property_base
    % Class of properties for linearly elastic transversely isotropic materials.
    %
    %
    
    properties
        % Young's moduli
        E1=[];      E2=[];  
        % Poisson ratio
        nu12=[];     nu23=[];
        % Shear moduli
        G13=[];     G23=[];
        % Thermal expansion coefficients
        alpha1= []; alpha2= []; alpha3= [];
        nconstrained =1;% how many tangent-moduli modes should be considered constrained?
    end
    
    properties (Hidden, SetAccess = private)
        % Young's moduli
           E3=[];
        % Poisson ratio
         nu13=[]; 
        % Shear moduli
        G12=[];  
        Precomputed_D = [];
        Precomputed_D_constrained=[];
        Precomputed_D_unconstrained=[];
    end
    
    methods
        
        function self = property_deformation_linear_transv_iso (Parameters)
            % Constructor.
            % Transversely isotropic material.
             % The local axis 1 is the axis of transverse isotropy.
            % Parameters:
            % required fields
            %     E1 =Young's modulus in the direction of the axis of transverse isotropy (the so-called longitudinal elastic modulus)
            %     E2 =Young's modulus in the plane of transverse isotropy (the so-called transverse elastic modulus)
            %     nu12 =Poisson ratio for loading along the axis of transverse isotropy
            %     nu23 =Poisson ratio for loading in the plane orthogonal to the axis of transverse isotropy
            %     G12=Shear modulus for sheer between the axis of transverse isotropy in the plane orthogonal to it
             %
             % The other parameters that define complete orthotropy are:
             % G23=E2/2/(1+nu23);
            % E3=E2;
            % nu13 = nu12;
            % G13=G12;
            %
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
            nu12=Parameters.nu12;
            E2=Parameters.E2;
            G12=Parameters.G12;
            nu23=Parameters.nu23;
            % The remaining coefficients are deduced using the assumption 
            % of the transverse isotropy
            G23=E2/2/(1+nu23);
            E3=E2;
            nu13 = nu12;
            G13=G12;
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
            % Prepare for the computation of the split material tangent moduli
            if ( isfield( Parameters, 'nconstrained'))
                self.nconstrained  = Parameters.nconstrained;
            end
            [DV,DL]=eig(self.Precomputed_D);
            [ignore, iix] =sort(diag(DL));;
            m1 = sqrt(3)*DV(:,iix(end-self.nconstrained+1:end));
            Iv =1/3*(m1 * m1');
            Id =eye(6)-Iv;
            self.Precomputed_D_constrained=Iv*self.Precomputed_D*Iv;
            self.Precomputed_D_unconstrained=Id*self.Precomputed_D*Id;
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
                    case 'constrained'
                        val = tangent_moduli_constrained(self, context);
                    case 'unconstrained'
                        val = tangent_moduli_unconstrained(self, context);
                end
            else
                val= self.Precomputed_D;
            end
            return;
        end
        
    end
    
    methods (Hidden, Access = private)
        
        function D = tangent_moduli_constrained(self, context)
            % Calculate the part of the material stiffness matrix that corresponds to
            % the single  constrained material direction.
            D= self.Precomputed_D_constrained;
        end
        
        function D = tangent_moduli_unconstrained(self, context)
            % Calculate the part of the material stiffness matrix that correspond to
            % the unconstrained material directions.
            D= self.Precomputed_D_unconstrained;
        end
        
        
    end
    
end
