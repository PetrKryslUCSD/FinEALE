classdef material_deformation_simo_j2_hard_triax < material_deformation_triax
    %     Class that represents deformable nonlinear J2 plasticity in the
    %     deformation gradient split format with mixed  kinematic/isotropic hardening.
    %
    % Reference
    % Interdisciplinary Applied Mathematics
    % Volume 7 1998
    % Computational Inelasticity
    % Authors: J. C. Simo, T. J. R. Hughes
    % ISBN: 978-0-387-97520-7 (Print) 978-0-387-22763-4 (Online)
    %
    
    properties (Access = private)
        mI = diag([1 1 1 0.5 0.5 0.5]);
        m1     = [1 1 1 0 0 0]';
        m1m1 =[1 1 1 0 0 0]'*[1 1 1 0 0 0]; %m1 * m1'
        bulk= 0;
        mu=0;
        eye3 =eye(3);
    end
    
    methods
        
        function self = material_deformation_simo_j2_hard_triax(Parameters)
            % Constructor.
            % Parameters:
            % none
            %
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if (nargin < 1)
                Parameters =struct( [] );
            end
            self = self@material_deformation_triax(Parameters);
            %             self.h = @(alpha)self.property.sigma_res - (self.property.sigma_res - self.property.sigma_y)*exp(-self.property.delta_expon*alpha) + self.property.K_lin*alpha;
            %             self.K = @(a)self.property.beta*self.h(a);
            %             self.H = @(a)(1-self.property.beta)*self.h(a);
            %             self.Dh = @(a)self.property.delta_expon*(self.property.sigma_res - self.property.sigma_y)*exp(-self.property.delta_expon*a) + self.property.K_lin;
            E = self.property.E;
            nu = self.property.nu;
            self.bulk = E  / 3 / (1 - 2*(nu));          % bulk modulus
            self.mu     = E / (2 * (1 + nu));           % shear modulus
        end
        
        function val = are_tangent_moduli_constant (self)
            % Is the material stiffness matrix independent of location (constant,
            % corresponding to a homogeneous material)?
            val =false;
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
            Fn1 = context.Fn1;
            Fn  = context.Fn;
            dt  = context.dt;
            
            % This is the case were the material state needs to be updated
            % and then the variable requested needs to be returned.
            
            sigma_y = self.property.sigma_y;
            sigma_res = self.property.sigma_res;
            delta_expon = self.property.delta_expon;
            K_lin = self.property.K_lin;
            beta = self.property.beta;
            SQRT_2_OVER_3= 0.81649658092772603272;
            
            % Incremental deformation gradient
            fn1=Fn1/Fn;
            fn1bar =((det(fn1))^(-1/3))*fn1;%
            bebarn1= fn1bar*ms.bebar*fn1bar';% trial Finger
            devbebarn1 =bebarn1-sum(diag(bebarn1))/3*self.eye3;
            % elastic trial
            sn1tr=self.mu*devbebarn1;
            strdevn = sqrt(sum(sum(sn1tr.*sn1tr)));% Double contraction
            Iebarn1 =sum(diag(bebarn1))/3;
            mubar=Iebarn1*self.mu;
            
            f1_trial=strdevn - SQRT_2_OVER_3 * K(ms.equiv_pl_def) ;
            if (f1_trial > 0)  % if plastic, correct
                n=sn1tr*(1/strdevn);    % unit normal
                %==================BOX3.1 Page 122 (Computational Inelasticity, Simo and Hughes) ===============
                Gamma = 0;   equiv_pl_def = ms.equiv_pl_def;
                Hprev =H(equiv_pl_def);
                Gamma = 0;  gv =g(equiv_pl_def);
                gv0=gv; gtol =max([gv0/1e6,100*eps]);
                iter=1;
                while abs(gv) > gtol %  Iteration loop
                    Dgv = Dg(equiv_pl_def + SQRT_2_OVER_3*Gamma) ;
                    Gamma = Gamma - gv/Dgv;
                    gv =g(equiv_pl_def + SQRT_2_OVER_3*Gamma) ;
                    iter=iter+1;
                end %  Iteration loop
                sn1 = sn1tr + (-2 * mubar * Gamma) * n;
                % update
                ms.equiv_pl_def = ms.equiv_pl_def + SQRT_2_OVER_3 * Gamma;
                ms.back_stress = ms.back_stress + SQRT_2_OVER_3*(H(ms.equiv_pl_def) - Hprev) * n;
            else
                sn1 = sn1tr;
            end
            % Update the intermediate configuration
            ms.bebar = sn1/self.mu + Iebarn1*self.eye3;
            % Compute the final stress
            Jn1 =  det(Fn1);
            pn1 = self.bulk/2*(Jn1-1/Jn1);
            sign1  =pn1*self.eye3 + sn1/Jn1;
            
            function value =h(alpha)
                value = sigma_res - (sigma_res - sigma_y)*exp(-delta_expon*alpha) + K_lin*alpha;
            end
            function value =K(alpha)
                value =beta*h(alpha);
            end
            function value =H(alpha)
                value =(1-beta)*h(alpha);
            end
            function value =Dh(alpha)
                value =delta_expon*(sigma_res - sigma_y)*exp(-delta_expon*alpha) + K_lin;
            end
            function value =g(alpha)
                value =-SQRT_2_OVER_3*K(alpha) + strdevn - (2*mubar*Gamma + SQRT_2_OVER_3*(H(alpha) - Hprev));
            end
            function value =Dg(alpha)
                Dhv =Dh(alpha);
                value= -2*mubar*(1 + (beta*Dhv+(1-beta)*Dhv)/3/mubar);
            end
            
            
            % This means to be programmed
            ms.deform_energy_density = 0;
            
            out =get_var();
            newms = ms;
            return;
            
            function out =get_var()
                switch context.output
                    case 'Cauchy'
                        out = self.stress_3x3t_to_6v(sign1);
                    case 'pressure'
                        out = -(sum(diag(sign1))/3);
                    case 'princCauchy'
                        [V,D]=eig(sign1);
                        out =sort(diag(D),'descend');
                    case {'vonMises','vonmises','von_mises','vm'}
                        stress = self.stress_3x3t_to_6v(sign1);
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
            Fn1 = context.Fn1;
            switch (stiff_type)
                case 'Eulerian'
                    ms=context.ms;
                    D=zeros(6,6);
                    dep=sqrt(eps);
                    ix=[1,4,5,;0,2,6;0,0,3];
                    ctx.Fn1 = context.Fn1; ctx.Fn = context.Fn;
                    ctx.dt = context.dt; ctx.dT = []; ctx.output='Cauchy';;
                    [sigmav0, ~] = state (self, ms, ctx);
                    Finc=zeros(3);
                    for i=1:3
                        for j=i:3
                            Finc(i,j)=Finc(i,j)+dep/2; Finc(j,i)=Finc(j,i)+dep/2;
                            ctx.Fn1 = Fn1 + Finc*Fn1;
                            [sigmav, ~] = state (self, ms, ctx);
                            D(:,ix(i,j)) =(sigmav-sigmav0)/(dep);
                            Finc(i,j)=0; Finc(j,i)=0;
                        end
                    end
                    D= (D+D')/2;
                    return;
                otherwise
                    error('Cannot handle');
                    D=[]; % return non-usable value
                    return;
            end
        end
        
        % Create a new material state.
        %
        function ms = newmatstate (self)
            ms.bebar  = eye(3);
            ms.back_stress = zeros(3);
            ms.equiv_pl_def = 0;
            ms.deform_energy_density = 0;
        end
        
    end
    
end
