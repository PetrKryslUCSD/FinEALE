classdef material_deformation_unrotated_j2 < material_deformation_triax
    %     Class that represents deformable nonlinear J2 plasticity in the
    %     unrotated format.
    %
    % Reference
    % WARP3D-­Release 17.6.0
    % 3-­D Dynamic Nonlinear Fracture Analyses of Solids  Using
    % Parallel Computers, Brian Healy, Arne Gullerud,  Kyle
    % Koppenhoefer, Arun Roy, Sushovan RoyChowdhury, Jason Petti, Matt 
    % Walters, Barron Bichon,  Kristine Cochran, Adam Carlyle,  James 
    % Sobotka,  Mark Messner  and Robert Dodds University	 of Illinois
    % at Urbana -­?Champaign, July 2015
    %
    
    properties
        D=[]; % The 6 x 6 constitutive matrix.
    end
    
    methods
        
        function self = material_deformation_unrotated_j2(Parameters)
            % Constructor.
            % Parameters:
            % none
            %
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if (nargin < 1)
                Parameters =struct( [] );
            end
            self = self@material_deformation_triax(Parameters);
            self.D = tangent_moduli(self.property, []);
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
            [R,U] = polardecomp(Fn1);
            
            if nargout == 1 % Pure retrieval,  no update
                out =get_var(); return
            end
            
            % This is the case were the material state needs to be updated
            % and then the variable requested needs to be returned.
            
            E = self.property.E;
            nu = self.property.nu;
            sigma_y = self.property.sigma_y;
            Hi = self.property.Hi;% Isotropic hardening modulus
            Hk = self.property.Hk;% Kinematic hardening modulus
            mu     = E / (2 * (1 + nu));                % shear modulus
            twomu = 2 * mu;
            SQRT_2_OVER_3= 0.81649658092772603272;
            lambda =E * nu /  (1 + nu) /  (1 - 2* (nu));
            dtl = lambda * dt; dttm = twomu * dt;
            
            midstep_F=compute_midstep_def_grad(Fn1,Fn);
            [midstep_R,~] = polardecomp(midstep_F);
            rel_F=compute_rel_def_grad(Fn1,Fn);
            [midstep_D,~]=compute_midstep_rate_of_def (rel_F, dt);
            
            t= ms.t; % previous unrotated stress
            prevt= t;
            d = midstep_R'*midstep_D*midstep_R;
            
            % compute the increment of unrotated stress
            lfact = lambda * dt;
            mfact = 2 * mu * dt;
            trd = sum(diag(d));
            incrt = mfact * d + lfact * trd*eye(3);
            
            % total unrotated stress
            t = t + incrt;
            
            % elastic predictor
            pressure = sum(diag(t))/ 3;
            tdev_predictor = t - pressure*eye(3);
            N= tdev_predictor;
            N = N - ms.back_stress;
            strdevn = sqrt(sum(sum(N.*N)));% Double contraction
            
            f1_trial=(strdevn - SQRT_2_OVER_3 * ms.sigma_y) ;
            if (f1_trial > 0)  % if plastic, correct
                N=N*(1/strdevn);    % normal unit
                % correct t deviator
                Gamma = (f1_trial / (2 * mu * (1 + (Hi+Hk) / 3 / mu)));
                t = tdev_predictor + (-2 * mu * Gamma) * N + pressure*eye(3);
                % update
                ms.equiv_pl_def = ms.equiv_pl_def + sqrt(2/3) * Gamma;
                ms.sigma_y = ms.sigma_y + sqrt(2/3) * Hi * Gamma;
                ms.back_stress = ms.back_stress + 2 / 3 * Hk * Gamma * N;
            end
            
            % save state: new unrotated stress
            ms.t = t;
            ms.deform_energy_density = ms.deform_energy_density ...
                + dt*0.5*(sum(sum(prevt.*midstep_D))+sum(sum(t.*midstep_D)))*det(Fn1);
            
            out =get_var();
            newms = ms;
            return;
            
           function out =get_var()
               % rotate to obtain the Cauchy stress
               Cauchystress= R*ms.t*R';
               switch context.output
                   case 'Cauchy'
                       out = self.stress_3x3t_to_6v(Cauchystress);
                   case 'pressure'
                       out = -(sum(diag(Cauchystress))/3);
                   case 'princCauchy'
                       [V,D]=eig(Cauchystress);
                       out =sort(diag(D),'descend');
                   case {'vonMises','vonmises','von_mises','vm'}
                       stress = self.stress_3x3t_to_6v(Cauchystress);
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
                    for i=1:3
                        for j=i:3
                            Finc=zeros(3); Finc(i,j)=Finc(i,j)+dep/2; Finc(j,i)=Finc(j,i)+dep/2;
                            ctx.Fn1 = Fn1 + Finc*Fn1;
                            [sigmav, ~] = state (self, ms, ctx);
                            D(:,ix(i,j)) =(sigmav-sigmav0)/dep;
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
            ms.sigma_y = self.property.sigma_y;
            ms.t = zeros(3);
            ms.back_stress = zeros(3);
            ms.equiv_pl_def = 0;
            ms.deform_energy_density = 0;
        end
        
    end
    
end
