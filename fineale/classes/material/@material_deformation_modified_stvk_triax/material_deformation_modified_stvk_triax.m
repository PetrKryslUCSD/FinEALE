classdef material_deformation_modified_stvk_triax < material_deformation_triax
% Class that represents deformable nonlinear hyperelastic 
% materials of the St.Venant-Kirchhoff type.
    %
    
    properties (Access = private)
        mI = diag([1 1 1 0.5 0.5 0.5]);
        m1     = [1 1 1 0 0 0]';
        m1m1 =[1 1 1 0 0 0]'*[1 1 1 0 0 0]; %m1 * m1'
        kappa= 0;
        lambda = 0;
        mu=0;
        eye3 =eye(3);
    end
    
    methods
        
        function self = material_deformation_modified_stvk_triax(Parameters)
            % Constructor.
            % Parameters:
            % none
            %
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if (nargin < 1)
                Parameters =struct( [] );
            end
            self = self@material_deformation_triax(Parameters);
            E = self.property.E;
            nu = self.property.nu;
            %             lambda = E * nu / (1 + nu) / (1 - 2*(nu));          % Lame constant #1
            self.mu     = E / (2 * (1 + nu));                        % shear modulus
            self.kappa  = E / 3 / (1 - 2*nu);                        % bulk modulus
            self.lambda = E * nu / (1 + nu) / (1 - 2*(nu));          % Lame constant #1
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
                        
            J=det(Fn1);
            % Cauchy-Green deformation tensor
            C=Fn1'*Fn1;
            %             S = self.kappa*log(J)*C^-1 + self.mu*(C- self.eye3);
            S =self.lambda*(J-1)*C^-1 + self.mu*(C- self.eye3);          % Lame constant #1
            
            out =get_var();
            newms = ms;
            return;
            
            function out =get_var()
                switch context.output
                    case 'Cauchy'
                        sigma=Fn1*(S/J)*Fn1';
                        out = self.stress_3x3t_to_6v(sigma);
                    case'strain_energy'
                        out = 1/2*S'*Egl;
                    case '2ndPK'
                        out= S;
                    case 'pressure'
                        sigma=Fn1*(self.stress_6v_to_3x3t(S)/J)*Fn1';
                        out = -(sum(diag(sigma))/3);
                    case {'vonMises','vonmises','von_mises','vm'}
                        sigma=Fn1*(self.stress_6v_to_3x3t(S)/J)*Fn1';
                        stress = self.stress_3x3t_to_6v(sigma);
                        s1=stress(1);s2=stress(2);s3=stress(3);
                        s4=stress(4);s5=stress(5);s6=stress(6);
                        out = sqrt(1/2*((s1-s2)^2+(s1-s3)^2+(s2-s3)^2+6*(s4^2+s5^2+s6^2)));
                    case 'equiv_pl_def'
                        out =ms.equiv_pl_def;
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
            switch (stiff_type)
                case 'Eulerian'
%                       ms=context.ms;
%                     D=zeros(6,6);
%                     dep=sqrt(eps);
%                     ix=[1,4,5,;0,2,6;0,0,3];
%                     ctx.Fn1 = context.Fn1; ctx.Fn = context.Fn;
%                     ctx.dt = context.dt; ctx.dT = []; ctx.output='Cauchy';;
%                     [sigmav0, ~] = state (self, ms, ctx);
%                     for i=1:3
%                         for j=i:3
%                             Finc=zeros(3); Finc(i,j)=Finc(i,j)+dep/2; Finc(j,i)=Finc(j,i)+dep/2;
%                             ctx.Fn1 = context.Fn1 + Finc*context.Fn1;
%                             [sigmav, ~] = state (self, ms, ctx);
%                             D(:,ix(i,j)) =(sigmav-sigmav0)/dep;
%                         end
%                     end
%                     D1= (D+D')/2;
                    %         Analytical expression:
                    C=context.Fn1'*context.Fn1;
                    Jn1= det(context.Fn1);
                    Cinv=inv(C);
                    ST  =self.lambda*(Jn1*tens4_ijkl(Cinv,Cinv) +(1-Jn1)*(tens4_ikjl(Cinv,Cinv)+tens4_iljk(Cinv,Cinv)))...
                        + self.mu*(tens4_identity()+tens4_transposor());;
                    D= tens4_symm_to_6(ST);
                    D = Lagrangean_to_Eulerian(self, D, context.Fn1);
                    
                    return;
                otherwise
                    error('Cannot handle');
                    D=[]; % return non-usable value
                    return;
            end
        end
        
        function val = are_tangent_moduli_constant (self)
        % Is the material stiffness matrix independent of location (constant, 
        % corresponding to a homogeneous material)?
            val =are_tangent_moduli_constant(self.property);
        end
             
    end
    
end

