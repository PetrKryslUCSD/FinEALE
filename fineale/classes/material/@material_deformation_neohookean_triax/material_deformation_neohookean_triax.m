classdef material_deformation_neohookean_triax < material_deformation_triax
    % Class that represents triaxial neohookean hyper elastic material.
    %
    
    properties
    end
    
    methods
        
        function self = material_deformation_neohookean_triax(Parameters)
            % Constructor.
            % Parameters:
            % none
            %
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if (nargin < 1)
                Parameters =struct( [] );
            end
            self = self@material_deformation_triax(Parameters);
        end
        
        function [out, newms] = update (self, ms, context)
            % Update material state.  Return the updated material state, and the
            % requested quantity (default is the stress).
            %   Call as:
            %     [out,newms] = update(m, ms, context)
            %  where
            %     m=material
            %     ms = material state
            %     context=structure
            %        with mandatory fields
            %           F= deformation gradient
            %        and optional fields
            %           output=type of quantity to output, and interpreted by the
            %           particular material; [] is returned when the material does not
            %           recognize the requested quantity to indicate uninitialized
            %           value.  It can be tested with isempty ().
            %              output ='Cauchy' - Cauchy stress; this is the default
            %                      when output type is not specified.
            %              output ='2ndPK' - 2nd Piola-Kirchhoff stress;
            %                  It is assumed that stress is output in 6-component vector form.
            %              output ='strain_energy'
            %   The output arguments are
            %     out=requested quantity
            %     newms=new material state; don't forget that if the update is final
            %           the material state newms must be assigned and stored.  Otherwise
            %           the material update is lost!
            %
            E = self.property.E;
            nu = self.property.nu;
            lambda = E * nu / (1 + nu) / (1 - 2*(nu));          % Lame constant #1
            mu     = E / (2 * (1 + nu));                        % shear modulus
            F1 = context.F;
            % Finger deformation tensor
            b=F1*F1';
            J=det(F1);
            sigma = (mu/J) * (b - eye(3,3)) + (lambda *log(J)/J) * eye(3,3);
            if isfield(context,'output')
                switch context.output
                    case 'Cauchy'
                        out = self.stress_3x3t_to_6v(sigma);
                    case'strain_energy'
                        C=F1'*F1;
                        out = mu/2*(trace(C)-3) - mu*log(J) + lambda/2*(log(J))^2;
                    case '2ndPK'
                        invF1=inv(F1);
                        S = J * (invF1*sigma*invF1');
                        out= stress_3x3t_to_6v(self.mater,S);
                    case 'pressure'
                        out = -(sum(diag(sigma))/3);
                    otherwise
                        out = stress_3x3t_to_6v(self.mater,sigma);
                end
            else
                out = stress;
            end
            newms = ms;
            return;
        end
        
        function D = tangent_moduli(self, context)
            % Calculate the material stiffness matrix.
            %   Call as:
            %     D = tangent_moduli(m, context)
            %   where
            %     m=material
            %    context=struct
            % with fields
            %  - mandatory
            %  - optional
            %     stiff_type=type of the stiffness,
            %        either 'Lagrangean'
            %        or 'Eulerian'
            %              requires field 'F'= current deformation gradient
            % the output arguments are
            %     D=matrix 6x6
            %
            stiff_type='Eulerian';
            if isfield(context,'stiff_type')
                stiff_type= context.stiff_type;
            end
            switch (stiff_type)
                case 'Eulerian'
                    E = self.property.E;
                    nu = self.property.nu;
                    J =det(context.F);
                    lambda = E * nu / (1 + nu) / (1 - 2*(nu));          % Lame constant #1
                    mu     = E / (2 * (1 + nu));                        % shear modulus
                    kappa  = E / 3 / (1 - 2*nu);                        % bulk modulus
                    mI = diag([1 1 1 0.5 0.5 0.5]);
                    m1     = [1 1 1 0 0 0]';
                    D = (lambda / J) * m1 * m1' + 2 * (mu - lambda*log(J))/J * mI;
                    return;
                otherwise
                    error('Cannot handle');
                    D=[]; % return non-usable value
                    return;
            end
        end
        
    end
    
end
