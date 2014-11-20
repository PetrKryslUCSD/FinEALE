% Class of elastic  – perfectly plastic material for triaxial strain/stress states.
%
% Describes material that is isotropically elastic / perfectly plastic.
% The Mises yield condition is used: hence, this is J-2 plasticity.
%
classdef material_deformation_ss_plastic_perf_j2_triax < material_deformation_linear_triax
    %
    
    properties
        % none
    end
    
    methods
        
        function self = material_deformation_ss_plastic_perf_j2_triax (Parameters)
            % Constructor.
            % Parameters:
            %     those recognized by material_deformation.
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if (nargin < 1)
                Parameters =struct( [] );
            end
            self = self@material_deformation_linear_triax(Parameters);
        end
        
        function D = tangent_moduli(self, context)
            % Calculate the material stiffness matrix.
            %
            % Arguments
            %     m=material
            %     context=structure with mandatory and optional fields that are required
            % by the specific material implementation.
            %
            % the output arguments are
            %     D=matrix 6x6 in the local material orientation Cartesian basis
            %
            %
            if (isfield(context,'strain'))
                Ev = context.strain;% strain in material coordinates
            else% This is an approximation valid only for small displacements
                gradu=context.F-eye(3);
                Ev = strain_3x3t_to_6v (self,(gradu+gradu')/2);
            end
            E = self.property.E;
            nu = self.property.nu;
            sigma_y = self.property.sigma_y;
            lambda = E * nu / (1 + nu) / (1 - 2*(nu));          % Lame constant #1
            mu     = E / (2 * (1 + nu));                        % shear modulus
            kappa  = E / 3 / (1 - 2*nu);                        % bulk modulus
            mI = diag([1 1 1 0.5 0.5 0.5]);
            m1     = [1 1 1 0 0 0]';
            eps1 = Ev;                        % total strain
            e1   = eps1 - 1/3*(eps1'*m1)*m1;             % deviatoric total strain
            s1_trial= stress_6v_to_3x3t (self,2 * mu * mI * (e1 - context.ms.eps_p));
            %     norm_s1_trial = sqrt(sum(s1_trial(1:3).^2 + 2*(s1_trial(4:6).^2)));
            norm_s1_trial = sqrt(trace(s1_trial*s1_trial));
            f1_trial = norm_s1_trial - sqrt(2/3)*sigma_y;  % trial yield f. value
            if (f1_trial > 0)
                n1 = s1_trial/norm_s1_trial;                        % unit normal to f1_trial
                % How come I'm not getting the same number?
                dgamma1 = f1_trial / (2 * mu);
                theta1 = 1 - 2*mu*dgamma1/norm_s1_trial;
                %         theta1 = 1 - f1_trial/norm_s1_trial;
                theta1bar =  theta1;
                n1=stress_3x3t_to_6v (self,n1);
                D = kappa*m1*m1' ...
                    + 2*mu*theta1*(mI-1/3*m1*m1')...
                    - 2*mu*theta1bar*n1*n1';
            else
                D = kappa*m1*m1' + 2*mu*(mI-1/3*m1*m1');
            end
        end
        
        function [out, newms] = update (self, ms, context)
            % Update material state.
            %
            % function [out, newms] = update (self, ms, context)
            %
            % Update material state.  Return the updated material state, and the
            % requested quantity (default is the stress).
            %   Call as:
            %     [out,newms] = update(m, ms, context)
            %  where
            %     m=material
            %     ms = material state
            %     context=structure
            %        with mandatory fields
            %           strain=strain vector  in the local material
            %               directions (which may be the same as the global coordinate
            %               directions)
            %        and optional fields
            %           output=type of quantity to output, and interpreted by the
            %               particular material; [] is returned when the material does not
            %               recognize the requested quantity to indicate uninitialized
            %               value.  It can be tested with isempty ().
            %                  output ='Cauchy' - Cauchy stress; this is the default
            %                      when output type is not specified.
            %                  output ='princCauchy' - principal Cauchy stress;
            %                  output ='pressure' - pressure;
            %                  output ='vol_strain' - volumetric strain;
            %                  output ='vonMises' - von Mises stress;
            %           outputRm=optional orientation matrix in which output should
            %               supplied
            %
            %   It is assumed that stress is output in 6-component vector form.
            %   The output arguments are
            %     out=requested quantity
            %           Remember: the output is expressed in the local material
            %           orientation  coordinates.
            %     newms=new material state; don't forget that if the update is final
            %           the material state newms must be assigned and stored.  Otherwise
            %           the material update is lost!
            %
            if (isfield(context,'strain'))
                Ev = context.strain;% strain in material coordinates
            else% This is an approximation valid only for small displacements
                gradu=context.F-eye(3);
                Ev = strain_3x3t_to_6v (self,(gradu+gradu')/2);
            end
            E = self.property.E;
            nu = self.property.nu;
            sigma_y = self.property.sigma_y;
            lambda = E * nu / (1 + nu) / (1 - 2*(nu));          % Lame constant #1
            mu     = E / (2 * (1 + nu));                        % shear modulus
            kappa  = E / 3 / (1 - 2*nu);                        % bulk modulus
            mI = diag([1 1 1 0.5 0.5 0.5]);
            m1     = [1 1 1 0 0 0]';
            eps1 = Ev;                        % total strain
            e1   = eps1 - 1/3*(eps1'*m1)*m1;             % deviatoric total strain
            s1_trial= stress_6v_to_3x3t (self,2 * mu * mI * (e1 - ms.eps_p));
            norm_s1_trial = sqrt(trace(s1_trial*s1_trial));
            f1_trial = norm_s1_trial - sqrt(2/3)*sigma_y;  % trial yield f. value
            sigma1 = kappa * (eps1'*m1) * eye(3,3) + s1_trial;       % trial stress
            if (f1_trial > 0)                                   % plastic correction
                n1 = s1_trial/norm_s1_trial;                    % unit normal to f1_trial
                dgamma1 = f1_trial / (2 * mu);                  %
                sigma1 = sigma1 - 2 * mu * dgamma1 * n1;        % correct stress
                ms.eps_p = ms.eps_p + dgamma1 * strain_3x3t_to_6v (self,n1);    % update material state
                ms.dgamma = dgamma1;
                ms.equiv_pl_def = ms.equiv_pl_def + sqrt(2/3) * ms.dgamma;
            end
            sigma1 =stress_3x3t_to_6v (self,sigma1);
            ms.strain_energy = ms.strain_energy + 1/2*sigma1'*(eps1 - ms.eps_p);
            stress = sigma1;
            context.ms =ms;
            tSigma = thermal_stress(self,context);
            stress = stress + tSigma;
            if isfield(context,'output')
                switch context.output
                    case 'Cauchy'
                        out = stress;
                    case 'vol_strain'
                        out = sum(Ev(1:3));
                    case 'pressure'
                        out = -(sum(stress(1:3))/3);
                    case 'princCauchy'
                        t = stress_6v_to_3x3t (self,stress);
                        [V,D]=eig(t);
                        out =sort(diag(D),'descend');
                    case {'vonMises','vonmises','von_mises','vm'}
                        s1=stress(1);s2=stress(2);s3=stress(3);
                        s4=stress(4);s5=stress(5);s6=stress(6);
                        out = sqrt(1/2*((s1-s2)^2+(s1-s3)^2+(s2-s3)^2+6*(s4^2+s5^2+s6^2)));
                    case 'equiv_pl_def'
                        out =ms.equiv_pl_def;
                    otherwise
                        out = [];
                end
            else
                out = stress;
            end
            newms = ms;;
        end
        
        % Create a new material state.
        %
        function v = newmatstate (self)
            v.eps_p = zeros(6,1); % plastic strains
            v.strain_energy = 0;  % energy of elastic deformation
            v.equiv_pl_def = 0; % equivalent plastic deformation
            v.dgamma = 0; % plastic multiplier (discrete version of the consistency parameter)
            return;
        end
        
    end
    
end

