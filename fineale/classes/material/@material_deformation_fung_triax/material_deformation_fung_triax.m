classdef material_deformation_fung_triax < material_deformation_triax
% Class that represents deformable nonlinear hyperelastic 
% materials of the Fung type.
    %
    
    properties
        c=1; % Constitutive parameter
        D=[]; % The 6 x 6 constant constitutive matrix.
    end
    
    properties (Access = private)
        Id=eye(3,3);
    end
    
    methods
        
        function self = material_deformation_fung_triax(Parameters)
            % Constructor.
            % Parameters:
            % none
            %
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if (nargin < 1)
                Parameters =struct( [] );
            end
            self = self@material_deformation_triax(Parameters);
            D = tangent_moduli(self.property, []);
            [V,L]=eig(D);
            self.c = max(diag(L));
            self.D=D/self.c;
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
            % Green-Lagrange strain
            Egl=self.strain_3x3t_to_6v(1/2*(C-self.Id));
            DEgl=self.D*Egl;
            Q =(DEgl'*Egl);
            S = (self.c*exp(Q))*DEgl;
            
            sigmav=self.stress_3x3t_to_6v(Fn1*(self.stress_6v_to_3x3t(S)/J)*Fn1');
            
            if isfield(context,'output')
                switch context.output
                    case 'Cauchy'
                        out = sigmav;
                    case'strain_energy'
                        out = 1/2*S'*Egl;
                    case '2ndPK'
                        out= S;
                    case 'pressure'
                        J=det(F1);
                        sigma=F1*(self.stress_6v_to_3x3t(S)/J)*F1';
                        out = -(sum(diag(sigma))/3);
                    otherwise
                        out = sigmav;
                end
            else
                out = S;
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
                    %                     ms=context.ms;
                    %                     D=zeros(6,6);
                    %                     dep=sqrt(eps);
                    %                     ix=[1,4,5,;0,2,6;0,0,3];
                    %                     ctx.Fn1 = context.Fn1; ctx.Fn = context.Fn;
                    %                     ctx.dt = context.dt; ctx.dT = []; ctx.output='Cauchy';;
                    %                     [sigmav0, ~] = state (self, ms, ctx);
                    %                     for i=1:3
                    %                         for j=i:3
                    %                             Finc=zeros(3); Finc(i,j)=Finc(i,j)+dep/2; Finc(j,i)=Finc(j,i)+dep/2;
                    %                             ctx.Fn1 = Fn1 + Finc*Fn1;
                    %                             [sigmav, ~] = state (self, ms, ctx);
                    %                             D(:,ix(i,j)) =(sigmav-sigmav0)/dep;
                    %                         end
                    %                     end
                    %                     D= (D+D')/2;
                    %
                    % Hyperelastic Modeling and Soft-Tissue Growth
                    % Integrated with the Smoothed Finite Element Method-SFEM,
                    % Von der Fakult¨at f¨ur Maschinenwesen
                    % der Rheinisch-Westf¨alischen Technischen Hochschule Aachen
                    % zur Erlangung des akademischen Grades eines Doktors der
                    % Ingenieurwissenschaften genehmigte Dissertation
                    % vorgelegt von
                    % Minh Tuan Duong,  2014
                    C=Fn1'*Fn1;
                    % Green-Lagrange strain
                    Egl=self.strain_3x3t_to_6v(1/2*(C-self.Id));
                    DEgl=self.D*Egl;
                    Q =(DEgl'*Egl);
                    D = (self.c*exp(Q))*self.D + (self.c*2*exp(Q))*DEgl*DEgl';
                    D = Lagrangean_to_Eulerian(self, D, Fn1);
                    return;
                otherwise
                    error('Cannot handle');
                    D=[]; % return non-usable value
                    return;
            end
        end
        
    end
    
end

