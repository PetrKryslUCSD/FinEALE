classdef property_base
% Base class for a group of classes that represent material properties.
%

    properties
        rho = [];%mass density
    end
    
    methods
    
        function self = property_base(Parameters)
        % Constructor.
        % Parameters:
        %     rho=mass density (default is 1.0)
        %
        % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
                   if (nargin < 1)
                return;
            end
            if isfield(Parameters,'rho')
                self.rho=Parameters.rho;
            end
        end
    
    end
    
end

