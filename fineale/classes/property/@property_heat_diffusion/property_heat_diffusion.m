classdef property_heat_diffusion < property_base
% Class for heat diffusion models of materials.
%

    properties 
        thermal_conductivity = [];% Thermal conductivity
        specific_heat= [];% Specific heat per unit volume
    end
    
    methods
    
        function self = property_heat_diffusion (Parameters)
 % Parameters:
%     thermal_conductivity=matrix of conductivities; 2x2 for 2-D domains, 3x3
%                  for 3-D domains; the components are given in the local
%                  coordinate frame (femmlock decides what that frame is).
%                  The conductivity may be supplied as a function handle, 
%                  that returns the conductivity is a function of temperature. 
%                  For instance kappa=@(T)(0.09+0.004*T)*eye(2);
%     specific_heat=specific heat per unit volume
%                  The specific heat coefficient may be supplied as 
%                  a function handle, that returns the conductivity 
%                  as a function of temperature. 
%                  For instance v=cp_function(T)
%     source=source of heat: function handle or a constant
%
% See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
        
           if (nargin < 1)
                Parameters =struct( [] );
            end
            self = self@property_base(Parameters);
            if isfield(Parameters,'thermal_conductivity')
                self.thermal_conductivity = Parameters.thermal_conductivity;
                if (~strcmp(class(self.thermal_conductivity),'function_handle'))
                    if (det(self.thermal_conductivity) < 0)
                        error('Negative determinant of thermal_conductivity!');
                    end
                end
            end
            if isfield(Parameters,'specific_heat')
                self.specific_heat=Parameters.specific_heat;
            end
        end
        
    end
    
end
