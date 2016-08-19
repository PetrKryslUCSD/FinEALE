classdef material_deformation_mooney_rivlin_triax < material_deformation_triax
% Class that represents deformable nonlinear hyperelastic 
% materials of the Fung type.
    %
    
    properties
        mu=0; % Constitutive parameter (equal to the shear modulus of linear elasticity)
        beta=0; % Constitutive parameter, zero means neohookean
        bulk=1; % Bulk modulus
    end
    
    properties (Access = private)
        mI = diag([1 1 1 0.5 0.5 0.5]);
        m1     = [1 1 1 0 0 0]';
        m1m1 =[1 1 1 0 0 0]'*[1 1 1 0 0 0]; %m1 * m1'
    end
    
    methods
        
        function self = material_deformation_mooney_rivlin_triax(Parameters)
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
            mu=E/2/(1+nu);
            self.mu = mu;
            if (isfield(Parameters, 'mu'))
                self.mu=Parameters.mu;
            end
            self.beta=0.0;
            if (isfield(Parameters, 'beta'))
                self.beta=Parameters.beta;
            end
            self.bulk =E/3/(1-2*nu);
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
            
            c1=1/2*self.mu*(1-self.beta);; c2=1/2*self.mu*(self.beta);
            
            % This formulation according to http://www.oofem.org/resources/doc/matlibmanual/html/node9.html
            b=Fn1*Fn1';% Cauchy-Green deformation tensor
            bdev =J^(-2/3)*b;%deviatoric part
            barI1b=trace(bdev);    barI2b=trace(bdev*bdev);
            sigma=(1/J)*(c1*(2*bdev-2/3*barI1b*eye(3)) +c2*(2*barI1b*b-4/3*barI2b*eye(3)-2*bdev*bdev)) + self.bulk *log(J)/J*eye(3,3);
            
            out =get_var();
            newms = ms;
            return;
            
            function out =get_var()
                switch context.output
                    case 'Cauchy'
                        out = self.stress_3x3t_to_6v(sigma);
                    case'strain_energy'
                        out = self.c01*(trace(b)-3) + self.c10*(trace(b^(-1)/det(b))-3) + (self.bulk *(J-1)^2)/2;
                    case '2ndPK'
                        invF1=inv(Fn1);
                        S = J * (invF1*sigma*invF1');
                        out= self.stress_3x3t_to_6v(S);
                    case 'pressure'
                        out = -(sum(diag(sigma))/3);
                    case 'princCauchy'
                        [V,D]=eig(sigma);
                        out =sort(diag(D),'descend');
                    case {'vonMises','vonmises','von_mises','vm'}
                        stress = self.stress_3x3t_to_6v(sigma);
                        s1=stress(1);s2=stress(2);s3=stress(3);
                        s4=stress(4);s5=stress(5);s6=stress(6);
                        out = sqrt(1/2*((s1-s2)^2+(s1-s3)^2+(s2-s3)^2+6*(s4^2+s5^2+s6^2)));
                    case 'equiv_pl_def'
                        out =0.0;
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
        
    end
    
end

