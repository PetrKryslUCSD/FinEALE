classdef material_deformation_ifr_j2 < material_deformation_triax
    % Class that represents deformable nonlinear material "instantaneous final
    % rotation" J2 plasticity.
    %
    % Reference
    %
    %
    
    properties
    end
    
    properties (Access = private)
        e=kinem_var_type_enum(); % Enumeration of kinematic quantities.
    end
    
    methods
        
        function self = material_deformation_ifr_j2(Parameters)
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
            sigma_y = self.property.sigma_y;
            Hi = self.property.Hi;% Isotropic hardening modulus
            Hk = self.property.Hk;% Kinematic hardening modulus
            mu     = E / (2 * (1 + nu));                % shear modulus
            bulk_mod  = E / 3 / (1 - 2*nu);             % bulk modulus
            twomu = 2 * mu;
            SQRT_2_OVER_3= 0.81649658092772603272;
            
            cauchy_prev = ms.cauchy;
            
            Fn1 = context.Fn1;
            Fn  = context.Fn;
            dt  = context.dt;  
            rel_F=compute_rel_def_grad(Fn1,Fn);
            [rel_R,rel_U] = polardecomp(rel_F);
            ln_rel_U=compute_ln_U (rel_U);
            
            rel_J = det(rel_F);
            ln_rel_J=log(rel_J);
            
            tr_ln_rel_U_3 =sum(diag(ln_rel_U)) / 3;
            ln_rel_U = ln_rel_U- tr_ln_rel_U_3*eye(3);
            
            %   /* elastic predictor: corotational */
            curr_press      = sum(diag(ms.cauchy))/3;
            press_predictor = curr_press + bulk_mod * ln_rel_J;
            dev_sig_corot= ms.cauchy - curr_press*eye(3) + twomu * ln_rel_U;
            
            alpha_corot = ms.back_stress;
            
            %   /* elastic predictor: after rotation */
            dev_sig_predictor=rel_R*dev_sig_corot*rel_R';
            alpha_predictor=rel_R*alpha_corot*rel_R';
            N=dev_sig_predictor - alpha_predictor;
            effective_stress = sum(sum(N.*N));% Double contraction
            effective_stress = sqrt(effective_stress);
            
            %   /* plasticity criterion */
            if (effective_stress > SQRT_2_OVER_3 * ms.curr_sigma_y)
                N= N / effective_stress;  %  /* normal unit */
                %     /* linear hardening */
                Gamma = (effective_stress - SQRT_2_OVER_3 * ms.curr_sigma_y) / (twomu * (1 + (Hi+Hk) / 3 / mu));
                dev_sig = dev_sig_predictor -twomu * Gamma*N;
                ms.equiv_pl_def = ms.equiv_pl_def + SQRT_2_OVER_3 * Gamma;
                alpha = alpha_predictor + N * (2 / 3 * Hk * Gamma);
                ms.curr_sigma_y = ms.curr_sigma_y+ SQRT_2_OVER_3 * Hi * Gamma;
            else
                %     /* elastic state */
                dev_sig= dev_sig_predictor;
                alpha= alpha_predictor;
            end
            
            %             /* save state */
            ms.cauchy = press_predictor*eye(3) + dev_sig;
            ms.back_stress= alpha;
            
            [msD,~]=compute_midstep_rate_of_def (rel_F, dt);
            ms.deform_energy_density = ms.deform_energy_density ...
                + dt*0.5*(sum(sum(cauchy_prev.*msD))+sum(sum(cauchy_prev.*ms.cauchy)))*det(Fn1);
            
            if isfield(context,'output')
                switch context.output
                    case 'Cauchy'
                        out = self.stress_3x3t_to_6v(ms.cauchy);
                    case 'pressure'
                        out = -(sum(diag(ms.cauchy))/3);
                    case 'princCauchy'
                        [V,D]=eig(ms.cauchy);
                        out =sort(diag(D),'descend');
                    case {'vonMises','vonmises','von_mises','vm'}
                        stress = self.stress_3x3t_to_6v(ms.cauchy);
                        s1=stress(1);s2=stress(2);s3=stress(3);
                        s4=stress(4);s5=stress(5);s6=stress(6);
                        out = sqrt(1/2*((s1-s2)^2+(s1-s3)^2+(s2-s3)^2+6*(s4^2+s5^2+s6^2)));
                    case 'equiv_pl_def'
                        out =ms.equiv_pl_def;
                    otherwise
                        out = [];
                end
            else
                out = self.stress_3x3t_to_6v(ms.cauchy);
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
            Fn1 = context.Fn1;
            switch (stiff_type)
                case 'Eulerian'
                    ms=context.ms;
                    D=zeros(6,6);
                    dep=sqrt(eps);
                    ix=[1,4,5,;0,2,6;0,0,3];
                    ctx.Fn1 = context.Fn1; ctx.Fn = context.Fn; 
                    ctx.dt = context.dt; ctx.dT = []; ctx.output='Cauchy';;
                    [sigmav0, ~] = update (self, ms, ctx);
                    for i=1:3
                        for j=i:3
                            Finc=zeros(3); Finc(i,j)=Finc(i,j)+dep/2; Finc(j,i)=Finc(j,i)+dep/2;
                            ctx.Fn1 = Fn1 + Finc*Fn1;
                            [sigmav, ~] = update (self, ms, ctx);
                            D(:,ix(i,j)) =(sigmav-sigmav0)/dep;
                        end
                    end
                    %                     for i=1:3
                    %                         for j=i:3
                    %                             Finc=zeros(3); Finc(i,j)=Finc(i,j)+dep/2; Finc(j,i)=Finc(j,i)+dep/2;
                    %                             ctx.Fn1 = Fn1 + Finc*Fn1;
                    %                             [sigmap, ~] = update (self, ms, ctx);
                    %                             ctx.Fn1 = Fn1 - Finc*Fn1;
                    %                             [sigman, ~] = update (self, ms, ctx);
                    %                             D(:,ix(i,j)) =(sigmap-sigman)/2/dep;
                    %                         end
                    %                     end
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
            ms.curr_sigma_y = self.property.sigma_y;
            ms.cauchy = zeros(3);
            ms.back_stress= zeros(3);
            ms.equiv_pl_def = 0;
            ms.deform_energy_density = 0;
        end
        
    end
    
end


