classdef material_deformation_linear_biax < material_deformation_biax
% Class of linear elastic material for biaxial strain/stress states.
%
% This class represents deformable materials, linear elastic,
% for biaxial strain/stress states.
%
%
    
    properties
     
    end
    
    methods
        
        function self = material_deformation_linear_biax (Parameters)
        % Constructor.
        % Parameters:
        %     those recognized by material_deformation, and
        % reduction = kind of model reduction
        % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if (nargin < 1)
                Parameters =struct( [] );
            end
            self = self@material_deformation_biax(Parameters);
            
        end
        
        function val = are_tangent_moduli_constant (self)
        % Is the material stiffness matrix independent of location (constant, 
        % corresponding to a homogeneous material)?
            val =are_tangent_moduli_constant(self.property);
        end
             
        function D = tangent_moduli(self, context)
        % Calculate the material stiffness matrix.
        %
        % function D = tangent_moduli(self, context)
        %
        %   Call as:
        %     D = tangent_moduli(m, context)
        %   where
        %     m=material
        %    context=structure with options
        %            interpreted by the property object
        %
        % the output arguments are
        %     D=matrix of size particular to the kind of model reduction in effect
        %
            D= self.property.tangent_moduli(context);
            switch self.reduction
                case 'axisymm'
                    D =D(1:4, 1:4);
                case 'strain'
                    D =D([1, 2, 4],[1, 2, 4]);
                case 'stress'
                    Dt =D(1:2, 1:2)-D(1:2,3)*D(3,1:2)/D(3,3);
                    D =D([1, 2, 4],[1, 2, 4]);
                    D(1:2, 1:2)= Dt;
                otherwise
                    error([' Reduction ' self.reduction ' not recognized']);
            end
        end
        
        function [out, newms] = update (self, ms, context)
        % Update material state.
        %
        % function [out, newms] = update (self, ms, context)
        %
        % Update material state.  Return the updated material state, and the
        % requested quantity (default is the stress).
        %   Call as:
        %     [out,newms] = update(m, ms, context)
        %  where
        %     m=material
        %     ms = material state
        %     context=structure
        %        with mandatory fields
        %           strain=strain vector  in the local material
        %               directions (which may be the same as the global coordinate
        %               directions)
        %        and optional fields
        %           output=type of quantity to output, and interpreted by the
        %               particular material; [] is returned when the material does not
        %               recognize the requested quantity to indicate uninitialized
        %               value.  It can be tested with isempty ().
        %                  output ='Cauchy' - Cauchy stress; this is the default
        %                      when output type is not specified.
        %                  output ='princCauchy' - principal Cauchy stress;
        %                  output ='pressure' - pressure;
        %                  output ='vol_strain' - volumetric strain;
        %           outputRm=optional orientation matrix in which output should 
        %               supplied   
        %
        %   It is assumed that stress is output in m-component vector 
        %           form. m=3 for plane stress, m=4 for plane strain or axially 
        %           symmetric.
        % The stress components are ordered in the "out" argument as:
        % 
        %   The output arguments are
        %     out=requested quantity
        %     newms=new material state; don't forget that if the update is final
        %           the material state newms must be assigned and stored.  Otherwise
        %           the material update is lost!
        %
            Ev = context.strain;
            D  = tangent_moduli (self, context);
            tSigma = thermal_stress(self,context);
            stress = D * Ev + tSigma;
            switch self.reduction
                case 'strain' % output through-the-thickness stress
                    alphas = self.property.thermal_expansion_coefficients ();
                    D6x6= self.property.tangent_moduli(context);
                    sz=D6x6(3,1:2)*Ev(1:2)-context.dT*D6x6(3,1:2)*alphas(1:2)-D6x6(3,3)*context.dT*alphas(3);
                    stress = [stress;sz];
            end
            if isfield(context,'output')
                switch context.output
                    case 'Cauchy'
                        out = stress;
                    case 'pressure'
                        switch self.reduction
                            case 'stress'
                                out = -sum(stress(1:2))/3;
                            case 'strain'
                                out = -sum(stress([1,2,4]))/3;
                            otherwise
                                out = -sum(stress(1:3))/3;
                        end
                    case 'vol_strain'
                        out = sum(Ev(1:2));%RAW this is incorrect for axi
                    case 'princCauchy'
                        switch self.reduction
                            case 'stress'
                                t = stress_3v_to_3x3t(self,stress);
                            case 'strain'
                                t = stress_4v_to_3x3t(self,stress);
                            otherwise %  axial symmetry
                                t = stress_4v_to_3x3t(self,stress([1,2,4,3]));
                        end
                        [V,D]=eig(t);
                        out =sort(diag(D),'descend');    
                    otherwise
                        out = [];
                end
            else
                out = stress;
            end
            newms = ms;
            return;
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
            alphas = self.property.thermal_expansion_coefficients ();
            switch self.reduction
                case 'axisymm'
                    D = tangent_moduli(self, context);
                    v = -D*context.dT*[alphas(1:3).*ones(3, 1); 0];
                case 'strain'
                    D=  self.property.tangent_moduli(context);% need 3-D
                    v = -context.dT*[D(1:2, 1:2)*(alphas(1:2).*ones(2,1))+...
                        alphas(3)*D(1:2,3); 0];
                otherwise
                    D = tangent_moduli(self, context);
                    v = -D*context.dT*[alphas(1:2).*ones(2, 1); 0];
            end
        end  
        
        function v = newmatstate (self)
        % Create a new material state.
        %
        % function v = newmatstate (self)
        %
        %   Call as:
        %     v = newmatstate(m)
        %   where
        %     m=an instance of a descendant of the class material
        %
            v= [];
        end
        
        function draw (self, gv, context)
        % Produce a graphic representation of data at the integration point.
        %
        % function draw (self, gv, context)
        %
        %
        % Input arguments
        % self = self
        % gv = graphic viewer
        % context = struct
        %  with mandatory fields
        %    xyz         - location at which the draw
        %    ms          - material state
        %    update_context     - context to be passed to the update () function.
        %  with optional fields
        %    quantity  - string: name of the quantity to display
        %                 if quantity=='stress'
        %                    the graphics representation is an ellipsoid of radii corresponding
        %                    to the principal values of stress;
        %                 otherwise
        %                    the graphics representation is a sphere of radius corresponding
        %                    to the value
        %    component - index that says which component of the `quantity' quantity;
        %                remember: the component points into the stress vector
        %    data_cmap - data color map to use to color the quantity
        %    es - transformation 3x3 matrix, with basis vectors in columns, that
        %         describes the local Cartesian system for the material point
        %
            quantity='stress';
            if (isfield(context,'quantity'))
                quantity= context.quantity;
            end
            scale = 1.0;
            if (isfield(context,'scale'))
                scale= context.scale;
            end
            switch (quantity)
                case 'stress'
                    [sigmav,ignoreme] = update(self, context.ms, context.update_context);
                    % sigmav=10*[10000; -100000; 200000];
                    % make sure I get a 3 x 3 tensor
                    if length (sigmav) ==3
                        sigmav (end+1) =0;
                    end
                    sigma=stress_4v_to_3x3t(self,sigmav);
                    if (~isempty(context.Rm))
                        if size(context.Rm)==[2 2]
                            Rm=[context.Rm zeros(2,1); zeros(1,2) 1];
                        else
                            Rm=context.Rm;
                        end
                        sigma=Rm*sigma*Rm';
                    end
                    [V,D]=eig(sigma);
                    component =1;
                    if isfield(context,'component')
                        component = context.component;
                    end
                    color =[1 1 0];
                    if isfield(context,'data_cmap')
                        color = map_data(context.data_cmap, sigmav(component));
                    end
                    context.facecolor= color;
                    draw_ellipsoid(gv, context.xyz, V, scale*diag(D), context);
                otherwise
                    %  do nothing;
            end
            return;
        end
        
    end
        
end
    
