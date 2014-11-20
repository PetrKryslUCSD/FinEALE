% Class of properties for linearly elastic isotropic materials. 
%
classdef property_deformation_linear_iso < property_base
    
    properties
        E = [];% Young's modulus
        nu = [];% Poisson ratio
        G = [];% shear modulus
        alpha = [];% coefficient of thermal expansion
    end
    
    properties (Hidden, SetAccess = private)
        Precomputed_D = [];
    end
    
    methods
        
        function self = property_deformation_linear_iso (Parameters)
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
            self = self@property_base(Parameters);
            if ( ~isfield( Parameters, 'E'))
                return
            end
            self.E = Parameters.E;
            self.nu =  0;
            self.G = [];% the default is to compute the shear modulus from the two above
            self.Precomputed_D =[];
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
            if isfield(Parameters,'G')
                self.G = Parameters.G;
                if (self.G < 0)
                    error('Negative shear modulus!');
                end
                self.nu = [];
            end
            self.alpha = 0.0;
            if isfield(Parameters,'alpha')
                self.alpha= Parameters.alpha;
            end
            % precompute
            E = self.E;
            nu = self.nu;
            if (isempty(nu))
                nu =E/self.G/2-1;
                if (self.nu < 0)
                    error('Negative Poisson ratio!');
                elseif (self.nu >= 0.5)
                    error('Incompressible material!');
                end
            end
            lambda = E * nu / (1 + nu) / (1 - 2*(nu));
            mu     = E / (2 * (1 + nu));
            mI = diag([1 1 1 0.5 0.5 0.5]);
            m1 = [1 1 1 0 0 0]';
            self.Precomputed_D = lambda * m1 * m1' + 2 * mu * mI;
        end
        
        function val = get.G (self)
        % Return the shear modulus.
            if (~ isempty(self.G))
                val = self.G;
            else
                val = self.E/2/(1+self.nu);
            end
        end
        
        function val = thermal_expansion_coefficients (self)
        % Return the coefficients of thermal expansion in all three directions.
            val = self.alpha*ones(3, 1);
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
        %                kind = string with possible values 'lambda',
        %                       'lambda_shear', 'bulk', 'bulk_shear'
        %
        % Output
        %     D=matrix 6x6 in the local material orientation Cartesian basis
        %
            if isfield(context,'kind')
                switch context.kind
                    case 'lambda'
                        val = tangent_moduli_lambda(self, context);
                    case 'lambda_shear'
                        val = tangent_moduli_lambda_shear(self, context);
                    case 'bulk'
                        val = tangent_moduli_bulk(self, context);
                    case 'bulk_shear'
                        val = tangent_moduli_bulk_shear(self, context);
                    case 'constrained'
                        %val = tangent_moduli_bulk(self, context);
                        val = tangent_moduli_lambda(self, context);
                    case 'unconstrained'
                        %val = tangent_moduli_bulk_shear(self, context);
                        val = tangent_moduli_lambda_shear(self, context);
                    otherwise
                        val = [];
                end
            else
                val= self.Precomputed_D;
            end
        end
        
    end
    
    
    methods (Hidden, Access = private)
        
        function D = tangent_moduli_lambda(self, context)
        % Calculate the part of the material stiffness matrix that corresponds to
        % the lambda Lame coefficient.
            E = self.E;
            nu = self.nu;
            lambda = E * nu / (1 + nu) / (1 - 2*(nu));
            m1 = [1 1 1 0 0 0]';
            D = lambda * m1 * m1';
            return;
        end
        
        function D = tangent_moduli_lambda_shear(self, context)
        % Calculate the part of the material stiffness matrix that correspond to shear.
        % Note: makes sense only for isotropic materials.
            E = self.E;
            nu = self.nu;
            mu     = E / (2 * (1 + nu));
            mI = diag([1 1 1 0.5 0.5 0.5]);
            D = 2 * mu * mI;
            return;
        end
        
        function D = tangent_moduli_bulk(self, context)
        % Calculate the part of the material stiffness matrix that corresponds to
        % the bulk modulus.
            E = self.E;
            nu = self.nu;
            B = E / 3 / (1 - 2*(nu));
            m1 = [1 1 1 0 0 0]';
            D = B * m1 * m1';
            return;
        end
        
        function D = tangent_moduli_bulk_shear(self, context)
        % Calculate the part of the material stiffness matrix that correspond to shear.
        % Note: makes sense only for isotropic materials.
            E = self.E;
            nu = self.nu;
            mu     = E / (2 * (1 + nu));
            D = mu * [2/3*[2 -1 -1; -1 2 -1; -1 -1 2] zeros(3,3); zeros(3,3) eye(3,3) ];
            return;
        end
        
    end
    
    % Some checking:
    %     mI = diag([1 1 1 0.5 0.5 0.5]);
    % m1 = [1 1 1 0 0 0]';
    % mId =eye(6) - 1/3*m1 * m1';
    % K  = sym('K','real');
    % G  = sym('G','real');
    % D = K*m1 * m1'+2*G*mI*mId
    %  (1/3*m1 * m1')*D*(1/3*m1 * m1')
    %  mId*D*mId
    %  (1/3*m1 * m1')*D*(1/3*m1 * m1')+ mId*D*mId-D
end
