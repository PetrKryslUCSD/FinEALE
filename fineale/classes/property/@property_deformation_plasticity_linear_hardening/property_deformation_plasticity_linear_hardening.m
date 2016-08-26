% Class of properties for linearly elastic isotropic materials. 
%
classdef property_deformation_plasticity_linear_hardening < property_deformation_linear_iso
    
    properties
        sigma_y =  []; % Yield stress
        Hi = 0.0;% isotropic hardening modulus
        Hk = 0.0;% kinematic hardening modulus
        Hn = 0.0;% hardening modulus
    end
    
    methods
        
        function self = property_deformation_plasticity_linear_hardening (Parameters)
            % Constructor.
            % Parameters:
            % Required fields
            %     E=Young's modulus
            %     nu=Poisson ratio
            % Optionally one may supply the shear modulus
            %     G=shear modulus
            %       The Poisson ratio will be then computed from Young's modulus and
            %       the shear modulus
            % Optional fields
            %     alpha = coefficient of thermal expansion, default 0.0
            %
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if (nargin < 1)
                Parameters =struct( [] );
            end
            self = self@property_deformation_linear_iso(Parameters);
            if ( ~isfield( Parameters, 'sigma_y'))
                return
            end
            self.sigma_y = 0.0;
            if isfield(Parameters,'sigma_y')
                self.sigma_y = Parameters.sigma_y;
                if (self.sigma_y <= 0)
                    error('Non-positive sigma_y!');
                end
            end
            if isfield(Parameters,'Hi')
                self.Hi = Parameters.Hi;
                if (self.Hi < 0)
                    error('Negative Hi!');
                end
            end
            if isfield(Parameters,'Hk')
                self.Hk = Parameters.Hk;
                if (self.Hk < 0)
                    error('Negative Hk!');
                end
            end
            if isfield(Parameters,'Hn')
                self.Hn = Parameters.Hn;
                if (self.Hn < 0)
                    error('Negative Hn!');
                end
            end
        end
        
    end

end
%
% function retobj = property_plastic_perf_j2_iso(varargin)
% 
% This class represents properties of perfect elastic-plastic materials
% whose behavior is described by J2  plasticity
% Always takes either zero or one argument.  For no arguments, the default
% object is created.  For one argument, the argument may be either an
% object of the same type, in which case the copy is created; or, it may be
% a struct with the following fields representing parameters.
% Parameters:
% required fields
%     as for property_linel_iso
%   where
%     E       - Young's modulus
%     nu      - Poisson ratio
%     sigma_y - yield stress
% optional fields
%     alpha = coefficient of thermal expansion
%
function retobj = property_plastic_j2_iso_DMPT2(varargin)
    class_name='property_plastic_j2_iso_DMPT2';
    if (nargin == 0)
        parent=property_linel_iso;
        self.sigma_y = 0.0;
        retobj = class(self,class_name, parent);
        return;
    elseif (nargin == 1)
        arg = varargin{1};
        if strcmp(class(arg),class_name) % copy constructor.
            retobj = arg;
            return;
        else
            Parameters =varargin{1};
            parent=property_linel_iso(Parameters);
            self.sigma_y = 0.0;
            self.Hi = 0.0;
            self.Hk = 0.0;
            self.Hn = 0.0;
            if isfield(Parameters,'sigma_y')
                self.sigma_y = Parameters.sigma_y;
                self.Hi = Parameters.Hi;
                self.Hk = Parameters.Hk;
                self.Hn = Parameters.Hn;
                if (self.sigma_y <= 0)
                    error('Non-positive sigma_y!');
                end
            end
            % precompute
             retobj = class(self,class_name, parent);
            return;
        end
    else
        error('Illegal arguments!');
    end
    return;
end
