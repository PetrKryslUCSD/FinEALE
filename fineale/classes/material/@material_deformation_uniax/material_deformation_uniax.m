% Class of material for uniaxial strain/stress states.
%
% This class represents deformable materials 
% for uniaxial strain/stress states.
%
classdef material_deformation_uniax < material_deformation
%
    
    properties
        % none
    end
    
    methods
        
        function self = material_deformation_uniax (Parameters)
        % Constructor.
        % Parameters:
        %     those recognized by material_deformation.
        % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if (nargin < 1)
                Parameters =struct( [] );
            end
            self = self@material_deformation(Parameters);
        end
        
        function cout=Lagrangean_to_Eulerian(self, C, F)
            % Convert a Lagrangean constitutive matrix to an Eulerian one.
            %
            % function cout=Lagrangean_to_Eulerian(self, C, F)
            %
            % Convert a Lagrangean constitutive matrix to an Eulerian one. 
            % NOTE: the Lagrangean matrix is presumed symmetric.
            % self = object
            % C    = Lagrangean constitutive matrix, 6x6, symmetric
            % F    = current deformation gradient, F_ij = \partial x_i / \partial X_j
            cout=C/det(F);
        end

        
        function T = stress_vector_rotation(self,Rm)
        % Calculate the rotation matrix for a stress vector.
        %
        %   T = stress_vector_rotation(self,Rm)
        %
        % Calculate the rotation of the stress vector to the
        % coordinate system given by the columns of the rotation matrix Rm.
        T =[Rm];
        end
        
        function Tbar = strain_vector_rotation(self,Rm)
        % Calculate the rotation matrix for a strain vector.
        %
        %   Tbar = strain_vector_rotation(self,Rm)
        %
        % Calculate the rotation of the strain vector to the
        % coordinate system given by the columns of the rotation matrix Rm.
        Tbar =[Rm];
        end
        
        function D = rotate_stiffness(self,D,Rm)
        % Rotate constitutive stiffness matrix of the material.
        %
        %         function D=transform_stiffness(self,D,Rm)
        %
        % Rotate constitutive stiffness matrix of the material to the
        % coordinate system given by the columns of the rotation matrix Rm.
        T =stress_vector_rotation(self,Rm);
        D = T*D*T';
        end
        
        function C = rotate_compliance(self,C,Rm)
        % Rotate constitutive compliance matrix of the material.
        %
        %   C = rotate_compliance(self,C,Rm)
        %
        % Rotate constitutive compliance matrix of the material to the
        % coordinate system given by the columns of the rotation matrix Rm.
        Tbar =strain_vector_rotation(self,Rm);
        C = Tbar*C*Tbar';
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
        
