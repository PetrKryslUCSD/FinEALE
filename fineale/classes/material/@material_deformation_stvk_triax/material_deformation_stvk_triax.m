classdef material_deformation_stvk_triax < material_deformation_triax
% Class that represents deformable nonlinear hyperelastic 
% materials of the St.Venant-Kirchhoff type.
    %
    
    properties
        D=[]; % The 6 x 6 constant constitutive matrix.
    end
    
    methods
        
        function self = material_deformation_stvk_triax(Parameters)
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
            Fn1 = context.Fn1;
            % Cauchy-Green deformation tensor
            C=Fn1'*Fn1;
            % Green-Lagrange strain
            Egl=self.strain_3x3t_to_6v(1/2*(C-eye(3,3)));;
            S = self.D*Egl;
            if isfield(context,'output')
                switch context.output
                    case 'Cauchy'
                        J=det(Fn1);
                        sigma=Fn1*(self.stress_6v_to_3x3t(S)/J)*Fn1';
                        out = self.stress_3x3t_to_6v(sigma);
                    case'strain_energy'
                        out = 1/2*S'*Egl;
                    case '2ndPK'
                        out= S;
                    case 'pressure'
                        J=det(Fn1);
                        sigma=Fn1*(self.stress_6v_to_3x3t(S)/J)*Fn1';
                        out = -(sum(diag(sigma))/3);
                    otherwise
                        J=det(Fn1);
                        sigma=Fn1*(self.stress_6v_to_3x3t(S)/J)*Fn1';
                        out = self.stress_3x3t_to_6v(sigma);
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
            switch (stiff_type)
                case 'Eulerian'
                    %                     ms=context.ms;
                    %                     D=zeros(6,6);
                    %                     dep=sqrt(eps);
                    %                     ix=[1,4,5,;0,2,6;0,0,3];
                    %                     ctx.F = context.F; ctx.dT = []; ctx.output='Cauchy';;
                    %                     [sigmav0, ~] = update (self, ms, ctx);
                    %                     for i=1:3
                    %                         for j=i:3
                    %                             Finc=zeros(3); Finc(i,j)=Finc(i,j)+dep/2; Finc(j,i)=Finc(j,i)+dep/2;
                    %                             ctx.F = context.F + Finc*(context.F);
                    %                             [sigmav, ~] = update (self, ms, ctx);
                    %                             D(:,ix(i,j)) =(sigmav-sigmav0)/dep;
                    %                         end
                    %                     end
                    %                     D= (D+D')/2;
                    D = Lagrangean_to_Eulerian(self, self.D, context.Fn1);
                    return;
                otherwise
                    error('Cannot handle');
                    D=[]; % return non-usable value
                    return;
            end
        end
        
    end
    
end

