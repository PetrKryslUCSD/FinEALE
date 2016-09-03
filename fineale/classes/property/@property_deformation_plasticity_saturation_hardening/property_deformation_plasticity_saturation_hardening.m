% Class of properties for elastoplastic isotropic materials with saturation
% nonlinear hardening.   The hardening is based on the equation
%
% h(a) =sigma_res - (sigma_res-sigma_y)*exp(-delta_expon*a) + K_lin*a
% K(a) = beta*h(a)
% H(a) = (1-beta)*h(a)
%
classdef property_deformation_plasticity_saturation_hardening < property_deformation_plasticity_linear_hardening
    
    properties
	  sigma_res = 0.0; % residual (asymptotic) yield stress
      delta_expon = 0.0;% 
      K_lin = 0.0;% linear isotropic hardening
      beta = 0.0;% hardening mixture parameter, beta = 0 means kinematic hardening, beta = 1 means isotropic hardening
    end
    
    methods
        
        function self = property_deformation_plasticity_saturation_hardening (Parameters)
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
            self = self@property_deformation_plasticity_linear_hardening(Parameters);
            if isfield(Parameters,'sigma_res')
                self.sigma_res = Parameters.sigma_res;
                if (self.sigma_res < 0)
                    error('Negative sigma_res!');
                end
            end
            if isfield(Parameters,'delta_expon')
                self.delta_expon = Parameters.delta_expon;
                if (self.delta_expon < 0)
                    error('Negative delta_expon!');
                end
            end
            if isfield(Parameters,'K_lin')
                self.K_lin = Parameters.K_lin;
                if (self.K_lin < 0)
                    error('Negative K_lin!');
                end
            end
            if isfield(Parameters,'beta')
                self.beta = Parameters.beta;
                if (self.beta < 0)
                    error('Negative beta!');
                end
                if (self.beta<0)|| (self.beta>1)
                error('Coefficient beta out of bounds <0,1>!');
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
            self.sigma_res = 0.0;
            self.delta_expon = 0.0;
            self.K_lin = 0.0;
            self.beta = 0.0;
            if isfield(Parameters,'sigma_y')
                self.sigma_y = Parameters.sigma_y;
                self.sigma_res = Parameters.sigma_res;
                self.delta_expon = Parameters.delta_expon;
                self.K_lin = Parameters.K_lin;
                self.beta = Parameters.beta;
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
