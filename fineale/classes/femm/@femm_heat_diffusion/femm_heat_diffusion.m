classdef femm_heat_diffusion < femm_base
    % Class for the finite element heat diffusion model.
    %
    %
    
    properties
        surface_transfer = 0;% surface heat transfer  coefficient
    end
    
    methods % constructor
        
        function self = femm_heat_diffusion (Parameters)
            % Constructor.
            % Parameters: those recognized by femm_base, plus
            %   surface_transfer = surface transfer coefficient;; optional, makes sense
            %           only for bounding surfaces.
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin <1
                Parameters = struct([]);
            end
            self = self@femm_base(Parameters);
            self.surface_transfer =0;
            if isfield(Parameters,'surface_transfer')
                self.surface_transfer= Parameters.surface_transfer;
            end
            switch self.fes.dim
                case 0
                case 1
                    if (isempty(self.Rm)),self.Rm=eye(1);end % for identity transformation
                case 2
                    if (isempty(self.Rm)),self.Rm=eye(2);end % for identity transformation
                case 3
                    if (isempty(self.Rm)),self.Rm=eye(3);end % for identity transformation
                otherwise
                    error (' Not implemented');
            end
        end
        
    end
    
end

