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
                gradu=context.Fn1-eye(3);
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
        
        function [out, newms] = state(self, ms, context)
            % Retrieve material state. 
            %
            %   function [out, newms] = state(self, ms, context)
            %
            % The method is used with either one output argument or with
            % two output arguments.
            % 1.  With only "out" as output argument:  Retrieve the
            %     the requested variable from the material state. 
            % 2.  With both output arguments: Update the material state,
            %     and return the requested variable from the material state 
            %     and the updated material state.
            %
            %     The requested quantity may or may not be supported by
            %     the particular material model.(default is the stress). 
            %
            %   Input arguments:
            % self=material
            % ms = material state
            % context=structure
            %    with mandatory fields
            %       Fn1= current deformation gradient (at time t_n+1)
            %       Fn=  previous converged deformation gradient (at time t_n)
            %    and optional fields
            %       output=type of quantity to output, and interpreted by the
            %           particular material; [] is returned when the material 
            %           does not recognize the requested quantity to indicate 
            %           uninitialized value.  It can be tested with isempty().
            %           output ='Cauchy' - Cauchy stress; this is the default
            %              when output type is not specified.
            %           output ='2ndPK' - 2nd Piola-Kirchhoff stress;
            %                  
            %              output ='strain_energy'
            %    It is assumed that stress is output in 6-component vector
            %    form. 
            %
            %   Output arguments:
            % out=requested quantity
            % newms=new material state; don't forget that if the update is
            %       final the material state newms must be assigned and
            %       stored.  Otherwise the material update is lost!
            %
            
            % Get the basic kinematic variables
            if (isfield(context,'strain'))
                Ev = context.strain;% strain in material coordinates
            else% This is an approximation valid only for small displacements
                gradu=context.Fn1-eye(3);
                Ev = strain_3x3t_to_6v (self,(gradu+gradu')/2);
            end
            
            if nargout == 1 % Pure retrieval,  no update
                out =get_var(); return
            end
            
            % This is the case were the material state needs to be updated
            % and then the variable requested needs to be returned.
            
            E = self.property.E;
            nu = self.property.nu;
            sigma_y = self.property.sigma_y;
            mu     = E / (2 * (1 + nu));                        % shear modulus
            kappa  = E / 3 / (1 - 2*nu);                        % bulk modulus
            mI = diag([1 1 1 0.5 0.5 0.5]);
            m1     = [1 1 1 0 0 0]';
            eps1 = Ev;                        % total strain
            e1   = eps1 - 1/3*(eps1'*m1)*m1;             % deviatoric total strain
            s1_trial= stress_6v_to_3x3t (self,2 * mu * mI * (e1 - ms.eps_p));
            norm_s1_trial = sqrt(trace(s1_trial*s1_trial)); % Effective stress
            f1_trial = norm_s1_trial - sqrt(2/3)*sigma_y;  % trial yield f. value
            sigma1 = kappa * (eps1'*m1) * eye(3,3) + s1_trial;       % trial stress
            if (f1_trial > 0)                                   % plastic correction
                n1 = s1_trial/norm_s1_trial;                    % unit normal to f1_trial
                dgamma1 = f1_trial / (2 * mu);                  %
                sigma1 = sigma1 - 2 * mu * dgamma1 * n1;        % correct stress
                ms.eps_p = ms.eps_p + dgamma1 * strain_3x3t_to_6v (self,n1);    % update material state
                ms.equiv_pl_def = ms.equiv_pl_def + sqrt(2/3) * dgamma1;
            end
            sigma1 =stress_3x3t_to_6v (self,sigma1);
            ms.deform_energy_density =  1/2*sigma1'*(eps1 - ms.eps_p);
            stress = sigma1;
            context.ms=ms;
            tSigma = thermal_stress(self,context);
            stress = stress + tSigma;
            ms.Cauchy=stress;
            
            out =get_var();
            newms = ms;
            return;
            
            function out =get_var()
                stress=ms.Cauchy;
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
                    case'strain_energy'
                        out = ms.deform_energy_density;
                    otherwise
                        out = [];
                end
            end
        end
        
        % Create a new material state.
        %
        function v = newmatstate (self)
            v.cauchy = zeros(6,1); % Cauchy stress
            v.eps_p = zeros(6,1); % plastic strains
            v.deform_energy_density = 0;  % energy of elastic deformation
            v.equiv_pl_def = 0; % equivalent plastic deformation
            return;
        end
        
    end
    
end

