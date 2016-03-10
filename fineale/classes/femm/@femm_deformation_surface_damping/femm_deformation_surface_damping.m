classdef femm_deformation_surface_damping < femm_deformation_linear
    % Class for the small displacement, small strain deformation
    % model for damping on surfaces exposed to fluid media.
    %
    % See discussion of constructors in <a href="matlab:helpwin 'FAESOR/classes/Contents'">classes/Contents</a>.
    % Options: those recognized by fem_deformation_linear plus
    %
    
    properties
        damping_abc_impedance=0;
    end
    
    methods % constructor
        
        function self = femm_deformation_surface_damping (Parameters)
            if nargin <1
                Parameters = struct([]);
            end
            self = self@femm_deformation_linear(Parameters);
            if (isfield(Parameters,'damping_abc_impedance'))
                self.damping_abc_impedance=Parameters.damping_abc_impedance;
            end
        end
        
    end
    
end
