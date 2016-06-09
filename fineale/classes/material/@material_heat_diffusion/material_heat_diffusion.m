classdef material_heat_diffusion < material_base
    % Class for heat diffusion models for materials.
    %
    % This class represents linear diffusion models in materials.
    %
    
    properties
        property= [];% heat diffusion property object
    end
    
    methods
        
        function self = material_heat_diffusion (Parameters)
            % Constructor.
            % Parameters:
            %     property=heat-diffusion property object
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin <1
                return
            end
            self.property = Parameters.property;
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
            %   with mandatory fields
            %     gradtheta= temperature gradient in the local material
            %       directions (which may be the same as the global coordinate
            %       directions)
            %     Rm = material orientation  matrix, may be supplied as empty  if it
            %       corresponds to an identity.
            %   and optional fields
            %     output=type of quantity to output, and interpreted by the
            %       particular material; [] is returned when the material does not
            %       recognize the requested quantity to indicate uninitialized
            %       value.  It can be tested with isempty ().
            %          output ='flux' - heat flux vector in the local material
            %                  directions; this is the default
            %                  when output type is not specified.
            %     outputRm = orientation matrix  of the coordinate system in which
            %       the output should be calculated;  if this matrix is not supplied,
            %       it is assumed that the output is to be provided  in the local
            %       material coordinate system; otherwise the output is first
            %       transformed to the global coordinate system, and then to the
            %       output coordinate system.
            %
            %   Output arguments:
            % out=requested quantity
            % newms=new material state; don't forget that if the update is
            %       final the material state newms must be assigned and
            %       stored.  Otherwise the material update is lost!
            %
           
            kappa = self.property.thermal_conductivity;
            gradtheta = reshape (context.gradtheta, [],1);% make column vector
            flux = - kappa * gradtheta;% in local material orientation coordinates
            if (isfield(context,'outputRm'))%  output coordinate system supplied?
                if (isfield(context,'Rm')) && ( ~isempty(context.Rm) )
                    flux = (context.Rm*flux);% in global coordinate system
                end
                flux = context.outputRm'*flux;% in output coordinate system
            end
            if isfield(context,'output')
                switch context.output
                    case 'flux'
                        out = flux;
                    otherwise
                        out = [];
                end
            else
                out = flux;
            end
            newms = ms;
            return;
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
            v = [];
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
            %    component - index that says which component of the `quantity' quantity;
            %                remember: the component points into the output vector
            %    data_cmap - data color map to use to color the quantity
            %    Rm - transformation 3x3 matrix, with basis vectors in columns, that
            %         describes the local Cartesian system for the material point
            %
            quantity='flux';
            if (isfield(context,'quantity'))
                quantity= context.quantity;
            end
            scale = 1.0;
            if (isfield(context,'scale'))
                scale= context.scale;
            end
            switch (quantity)
                case 'flux'
                    flux = state(self, context.ms, context.update_context);
                    if ( isfield(  context,'Rm' ))
                        if (~isempty(context.Rm))
                            flux=context.Rm*flux;
                        end
                    end
                    component =1;
                    if isfield(context,'component')
                        component = context.component;
                    end
                    color =[1 1 0];
                    if isfield(context,'data_cmap')
                        color = map_data(context.data_cmap, flux(component));
                    end
                    context.facecolor= color;
                    draw_arrow(gv, context.xyz, scale*flux, context);
                otherwise
                    %  do nothing;
            end
            return;
        end
        
        
    end
    
end
