% Class of linear elastic material for triaxial strain/stress states.
%
% This class represents deformable materials, linear elastic,
% for triaxial strain/stress states.
%
classdef material_deformation_linear_triax < material_deformation_triax
%
    
    properties
        % none
    end
    
    methods
        
        function self = material_deformation_linear_triax (Parameters)
        % Constructor.
        % Parameters:
        %     those recognized by material_deformation.
        % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if (nargin < 1)
                Parameters =struct( [] );
            end
            self = self@material_deformation_triax(Parameters);
        end
        
        function val = are_tangent_moduli_constant (self)
        % Is the material stiffness matrix independent of location (constant, 
        % corresponding to a homogeneous material)?
            val =are_tangent_moduli_constant(self.property);
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
            D= self.property.tangent_moduli(context);
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
            
            % This is the case were the material state needs to be updated
            % and then the variable requested needs to be returned.
            
            D  = tangent_moduli (self, context);% local material stiffness
            tSigma = thermal_stress(self,context);% stress in local coordinates
            stress = D * Ev + tSigma;
            
            out =get_var();
            newms = ms;
            return;
            
            function out =get_var()
                switch context.output
                    case 'Cauchy'
                        out = stress;
                    case 'pressure'
                        out = -(sum(stress(1:3))/3);
                    case 'princCauchy'
                        [V,D]=eig(self.stress_6v_to_3x3t(stress));
                        out =sort(diag(D),'descend');
                    case {'vonMises','vonmises','von_mises','vm'}
                        s1=stress(1);s2=stress(2);s3=stress(3);
                        s4=stress(4);s5=stress(5);s6=stress(6);
                        out = sqrt(1/2*((s1-s2)^2+(s1-s3)^2+(s2-s3)^2+6*(s4^2+s5^2+s6^2)));
                    case 'vol_strain'
                        out = sum(Ev(1:3));
                    case 'equiv_pl_def'
                        out =0.0;
                    case'strain_energy'
                        out = (stress-tSigma)'*Ev/2;
                    otherwise
                        out = [];
                end
            end
        end
        
        function v = thermal_stress(self,context)
        % Calculate vector of thermal stress components.
        %
        % function v = thermal_stress(self,context)
        %
        %   Call as:
        %     v = thermal_stress(m,context)
        %  where
        %     m=material
        %     context=structure; see the update() method
        %
            if ( ~isempty(context.dT))
                D = tangent_moduli(self, context);
                alphas = self.property.thermal_expansion_coefficients();
                v =-D*context.dT*[alphas.*ones(3, 1); zeros(3, 1)];
            else
                v=zeros(6, 1);
            end
        end
        
    end
        
end
    
    
    
        %         function draw (self, gv, context)
        %         % Produce a graphic representation of data at the integration point.
        %         %
        %         % function draw (self, gv, context)
        %         %
        %         %
        %         % Input arguments
        %         % self = self
        %         % gv = graphic viewer
        %         % context = struct
        %         %  with mandatory fields
        %         %    xyz         - location at which the draw
        %         %    ms          - material state
        %         %    update_context     - context to be passed to the update () function.
        %         %  with optional fields
        %         %    quantity  - string: name of the quantity to display
        %         %                 if quantity=='stress'
        %         %                    the graphics representation is an ellipsoid of radii corresponding
        %         %                    to the principal values of stress;
        %         %                 otherwise
        %         %                    the graphics representation is a sphere of radius corresponding
        %         %                    to the value
        %         %    component - index that says which component of the `quantity' quantity;
        %         %                remember: the component points into the stress vector
        %         %    data_cmap - data color map to use to color the quantity
        %         %    es - transformation 3x3 matrix, with basis vectors in columns, that
        %         %         describes the local Cartesian system for the material point
        %         %
        %             quantity='stress';
        %             if (isfield(context,'quantity'))
        %                 quantity= context.quantity;
        %             end
        %             scale = 1.0;
        %             if (isfield(context,'scale'))
        %                 scale= context.scale;
        %             end
        %             switch (quantity)
        %                 case 'stress'
        %                     [sigmav,ignoreme] = update(self, context.ms, context.update_context);
        %                     [V,D]=eig(stress_6v_to_3x3t(self,sigmav));
        %                     component =1;
        %                     if isfield(context,'component')
        %                         component = context.component;
        %                     end
        %                     color =[1 1 0];
        %                     if isfield(context,'data_cmap')
        %                         color = map_data(context.data_cmap, sigmav(component));
        %                     end
        %                     context.facecolor= color;
        %                     draw_ellipsoid(gv, context.xyz, V, scale*diag(D), context);
        %                 otherwise
        %                     %  do nothing;
        %             end
        %             return;
        %         end
        
