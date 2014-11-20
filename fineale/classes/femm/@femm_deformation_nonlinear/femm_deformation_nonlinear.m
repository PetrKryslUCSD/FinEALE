classdef femm_deformation_nonlinear < femm_deformation_linear
    %    Class for large displacement, large strain mechanics to be used with
    %     continuum, isoparametric elements.
    
    properties
        matstates= {};% Material states, for each element and each quadrature point
    end
    
    
    methods % constructor
        
        function self = femm_deformation_nonlinear (Parameters)
            % Constructor.
            % Parameters: those recognized by model_base plus
            %    N = Outer normal, either as an array of numbers,or as a function handle.
            %     Signature as in this example (tangent is in matrix with tangent
            %     vectors as columns).
            %     function n= ABC_normal(x,tangent)
            %         n= cross(tangent(:,1),tangent(:,2));
            %         n=n/norm(n);
            %     end
            %
            % See discussion of constructors in <a href="matlab:helpwin 'fineale/classes/Contents'">classes/Contents</a>.
            if nargin <1
                Parameters = struct([]);
            end
            self = self@femm_deformation_linear(Parameters);
            self.matstates = {};
            if isempty(self.fes)||(count(self.fes) < 1)
                error('Need at least one finite element!');
            end
            integration_rule = self.integration_rule;
            if ~isempty(integration_rule)
                npts_per_fe = integration_rule.npts;
                self.matstates = cell(count(self.fes),npts_per_fe);
                for i=1:count(self.fes)
                    for j=1:npts_per_fe;
                        self.matstates{i,j} = newmatstate(self.material);
                    end
                end
            else
                self.matstates = {};
            end
        end
        
    end
    
    methods
        
        
    end
    
    
end


% function retobj = feblock_defor_nonlinear (varargin)
% class_name='feblock_defor_nonlinear';
% if (nargin == 0)
% parent = feblock;
% self.hBlmat= [];
% self.matstates = [];
% self.Rm= [];
% retobj = class(self,class_name,parent);
% return;
% elseif (nargin == 1)
% arg = varargin{1};
% if strcmp(class(arg),class_name) % copy constructor.
% retobj = arg;
% return;
% else
% options =varargin{1};
% parent = feblock(varargin{:});
% self.hBlmat= [];
% self.matstates = [];
% gcells = get(parent,'gcells');
% if (count(gcells) < 1)
% error('Need at least one gcell!');
% end
% integration_rule = get(parent,'integration_rule');
% if ~isempty(integration_rule)
% npts_per_gcell = get(integration_rule,'npts_per_gcell');
% self.matstates = cell(count(gcells),npts_per_gcell);
% for i=1:count(gcells)
% for j=1:npts_per_gcell;
% self.matstates{i,j} = newmatstate(get(parent,'mater'));
% end
% end
% else
% self.matstates = [];
% end
% dim= get(gcells,'dim');
% switch dim
% case 1
% self.hBlmat=@feblock_defor_nonlinear_Blmat1;
% case 2
% if get(gcells,'axisymm')
% self.hBlmat=@feblock_defor_nonlinear_Blmat2axisymm;
% else
% self.hBlmat=@feblock_defor_nonlinear_Blmat2;
% end
% case 3
% self.hBlmat=@feblock_defor_nonlinear_Blmat3;
% otherwise
% error (' Not implemented');
% end
% retobj = class(self,class_name,parent);
% return;
% end
% else
% error('Illegal arguments!');
% end
% return;
% end
