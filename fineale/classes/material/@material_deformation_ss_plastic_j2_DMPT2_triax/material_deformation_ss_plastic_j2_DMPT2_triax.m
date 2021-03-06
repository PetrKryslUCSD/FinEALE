% Class of elastic  � perfectly plastic material for triaxial strain/stress states.
%
% Describes material that is isotropically elastic / perfectly plastic.
% The Mises yield condition is used: hence, this is J-2 plasticity.
%
classdef material_deformation_ss_plastic_j2_DMPT2_triax < material_deformation_linear_triax
    %
    
    properties
        % none
    end
    
    methods
        
        function self = material_deformation_ss_plastic_j2_DMPT2_triax (Parameters)
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
            Hi = self.property.Hi;
            Hk = self.property.Hk;
            mu     = E / (2 * (1 + nu));                        % shear modulus
            kappa  = E / 3 / (1 - 2*nu);                        % bulk modulus
            eps1 = Ev;                        % total strain
            % Outside call
            [ignore, ignore, D ] = DMPT2(self,  Hi, Hk, sigma_y, eps1, context.ms, 2*mu, kappa);
            
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
            Hi = self.property.Hi;
            Hk = self.property.Hk;
            %
            EE  = E;
            ni  = nu;
            sy0 = sigma_y;
            G     = EE/(2*(1+ni));        % Elastic material shear modulus
            KK    = EE/(3*(1-2*ni));      % Elastic material bulk modulus
            twoG = 2*G;
            eps1 = Ev;                        % total strain
            % chiamata
            [ms, sigma1, D ] = DMPT2(self,  Hi, Hk, sy0, eps1, ms, 2*G, KK);
            ms.strain_energy =  1/2*sigma1'*(eps1 - ms.eps_p);
            stress = sigma1;
            context.ms =ms;
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
            v.Cauchy = zeros(6,1); % Cauchy stress
            v.eps_p = zeros(6,1); % plastic strains
            v.eps_p9 = zeros(9,1); % plastic strains
            v.strain_energy = 0;  % energy of elastic deformation
            v.equiv_pl_def = 0; % equivalent plastic deformation
            v.dgamma = 0; % plastic multiplier (discrete version of the consistency parameter)
            v.alpha = zeros(6,1); % backstress
            v.eps_n = zeros(6,1); % total strain at previously converged step t_n
            return;
        end
        
    end
    
    methods (Access = private)
        
        % function [sig,ep,alp,gamma,Cdisc,lam2]=DMPT2(eps_n,eps,eep_n,alfa,gamman,E)
        function [ms, sig6, Cdisc ] = DMPT2(self, Hi, Hk, sy0, eps, ms, twoG, KK)
            %
            M1 = eye(6);
            M1(4:6,4:6) = 1/2*M1(4:6,4:6);
            eps_n = M1*ms.eps_n; % total strain at t_n
            eep_n = M1*ms.eps_p;
            alfa = ms.alpha;
            gamman = ms.dgamma;
            % G     = EE/(2*(1+ni));   % Elastic material shear modulus
            % twoG  = 2*G;
            % KK    = EE/(3*(1-2*ni)); % Elastic material bulk modulus
            % Material parameters
            tol=1.0e-12;
            k1    = twoG/(twoG+Hk+Hi);
            k2    =    1/(twoG+Hk+Hi);
            al    = 1/2;
            k3    = 1/al;
            lam1=0;
            lam2=0;
            %
            % Compute deviatoric quantities @ t_n+1
            %
            eps = M1*eps;
            theta       =     eps(1,1)+eps(2,1)+eps(3,1); % volumetric strain
            epsd(1,1)   =     eps(1,1) - theta/3;         % e11 dev. component
            epsd(2,1)   =     eps(2,1) - theta/3;         % e22 dev. component
            epsd(3,1)   =     eps(3,1) - theta/3;         % e33 dev. component
            epsd(4:6,1) =     eps(4:6,1);    % [e12,e23,e31] dev. components
            %
            % @ t_n
            %
            thetan       =     eps_n(1,1) + eps_n(2,1) + eps_n(3,1); % volumetric strain
            epsnd(1,1)   =     eps_n(1,1) - thetan/3;         % e11 dev. component
            epsnd(2,1)   =     eps_n(2,1) - thetan/3;         % e22 dev. component
            epsnd(3,1)   =     eps_n(3,1) - thetan/3;         % e33 dev. component
            epsnd(4:6,1) =     eps_n(4:6,1);
            %
            % sigma deviatoric @ t_n
            %
            ssn = twoG*(epsnd - eep_n);
            %
            UNI  = [1 1 1 0 0 0];         % 2nd-order identity tensor
            Iden = eye(6);                % 4th-order identity tensor
            % Idev = Iden - 1/3*(UNI'*UNI); % 4th-order deviatoric tensor
            Idev = M1 - 1/3*(UNI'*UNI); % 4th-order deviatoric tensor
            M2 = Iden;
            M2(4:6,4:6) = 2*M2(4:6,4:6);
            %
            % elastic trial @ t_n+al
            %
            epsdal=al*epsd+(1-al)*epsnd; % deviatoric strain @ t_n+al
            epal=eep_n;
            ssal=twoG*(epsdal-eep_n);
            pp = KK*(al*theta+(1-al)*thetan); % elastic pressure at t_n+al
            alpal=alfa;
            Sigaltr=ssal-alpal;
            Sigaltr_norm = normV(Sigaltr);
            Sigal=Sigaltr;
            gammaal=gamman;
            nnal=Sigaltr/Sigaltr_norm;
            %
            nnal_diad_nnal = (nnal*nnal');
            %     nnal_diad_nnal(:,4:6) = 2*nnal_diad_nnal(:,4:6);
            %             nnal_diad_nnal(:,4:6) = nnal_diad_nnal(:,4:6);
            %
            F=sy0+Hi*gammaal;
            dif=Sigaltr_norm-F;
            Ddisc_1 = zeros(6,6);
            g11 = 0;
            g12 = 0;
            if dif/F > tol;
                % STEP 1: plastic
                %     disp('step1..plastic')
                lam1=dif/(twoG+Hi+Hk);
                Ylam1=(twoG+Hk)*lam1;
                % update material state @ t_n+al
                epal=epal+lam1*nnal;
                alpal=alpal+Hk*lam1*nnal;
                ssal=ssal-twoG*lam1*nnal;
                Sigal=Sigaltr-Ylam1*nnal;
                gammaal=gammaal+lam1;
                % compute Ddisc_1
                g11     = (twoG*al*lam1/Sigaltr_norm);
                g12     = (k1*al - twoG*al*lam1/Sigaltr_norm);
                %
                Ddisc_1 = g11 * Idev + g12 * nnal_diad_nnal;
            end
            %
            ms.dgamma = gammaal;
            ms.eps_p9 = strain_6v_to_9v(self,eep_n);
            ms.alpha = alfa;
            ms.eps_p = M2*eep_n; %strain_9v_to_6v(self,ei);
            ms.equiv_pl_def = sqrt(2/3) * ms.dgamma;
            ms.eps_n = al*M2*eps+(1-al)*M2*eps_n;
            sig6(1:3,1) = ssal(1:3,1) + pp; % total stress at t_n+1
            sig6(4:6,1) = ssal(4:6,1);      % total stress at t_n+1
            %
            % check consistency at t_n+a
            % disp(['consistency check at t_n+a: ',num2str((normV(Sigal)-(sy0+gammaal*Hi))/(sy0+gammaal*Hi))])
            % STEP 2
            Sign = ssn-alfa;
            ep=(1/al)*epal-((1-al)/al)*eep_n;
            ss=ssal/al - ((1-al)/al)*ssn;
            alp=alpal/al- ((1-al)/al)*alfa;
            Sigtr=Sigal/al - ((1-al)/al)*Sign;
            Sigtr_norm = normV(Sigtr);
            gamma=gammaal/al - ((1-al)/al)*gamman;
            nn=Sigtr/Sigtr_norm;
            %
            nn_diad_nn = (nn*nn');
            % nn_diad_nn(:,4:6) = 2*nn_diad_nn(:,4:6);
            nn_diad_nn(:,4:6) = nn_diad_nn(:,4:6);
            
            nn_diad_nnal = (nn*nnal');
            % nn_diad_nnal(:,4:6) = 2*nn_diad_nnal(:,4:6);
            nn_diad_nnal(:,4:6) = nn_diad_nnal(:,4:6);
            %
            Ddisc_2 = zeros(6,6);
            if (Sigtr_norm-(sy0+gamma*Hi))/(sy0+gamma*Hi)>tol;
                % STEP 2: plastic
                %     disp('step2..plastic')
                lam2=(Sigtr_norm-(sy0+gamma*Hi))/(twoG+Hk+Hi);
                % update material state @ t_n+1
                ep=ep+lam2*nn;
                alp=alp+Hk*lam2*nn;
                ss=ss-twoG*lam2*nn;
                Sig = ss-alp;
                gamma=gamma+lam2;
                % consistency check at t_n+1
                %     disp(['consistency check at t_n+1: ',num2str((normV(Sig)-(sy0+gamma*Hi))/(sy0+gamma*Hi))])
                % compute Ddisc_2
                a1 = k1*k2*Hi;
                a2 = k2*k3*(twoG*al -(twoG+Hk)*g11);
                a3 = k2*k3*(twoG+Hk)*g12*(nn'*M2*nnal);
                %
                a4 = k3/Sigtr_norm*(twoG*al -(twoG+Hk)*g11)*lam2;
                a5 = k3/Sigtr_norm*(twoG+Hk)*(nn'*M2*nnal)*g12*lam2;
                a6 = a4;
                a7 = k3/Sigtr_norm*(twoG+Hk)*g12*lam2;
                %
                g21 = a4;
                g22 = a7;
                g23 = -(a1+a3+a5);
                g24 = a2-a4;
                %
                Ddisc_2 = g21 * Idev         + g22 * nnal_diad_nnal + ...
                    g23 * nn_diad_nnal + g24 * nn_diad_nn ;
            end
            ms.dgamma = gamma;
            ms.eps_p9 = strain_6v_to_9v(self,ep);
            ms.alpha = alp;
            ms.eps_p = M2*ep;
            ms.equiv_pl_def = sqrt(2/3) * ms.dgamma;
            % sig6 = sig;%stress_9v_to_6v(self,sig);
            pp = KK*theta; % elastic pressure at t_n+1
            sig6(1:3,1) = ss(1:3,1) + pp; % total stress at t_n+1
            sig6(4:6,1) = ss(4:6,1);      % total stress at t_n+1
            % consistent tangent moduli
            dsde = twoG*Idev - twoG/al * Ddisc_1 - twoG * Ddisc_2;
            Cdisc = dsde + KK*(UNI'*UNI);
            
            return
            
            function n=normV(w)% deve essere un vettore colonna
                n = 0;
                %----------------------------
                for i = 1:3
                    n = n + w(i,1)^2;
                    n = n + 2*w(i+3,1)^2;
                end
                n = sqrt(n);
            end
        end

    end
    
    
end
