% Get property from the specified object and return the value.
%
% function val = get(self,varargin)
%
% See discussion of the get/set methods in <a href="matlab:helpwin
% 'sofea/classes/Contents'">classes/Contents</a>. 
%
function val = get(self,varargin)
if (nargin == 1)
    val = {{{'npts'}, {'number of quadrature points, scalar'}};
        {{'npts_per_gcell'}, {'total number of quadrature points per gcell (same as npts), scalar'}};
        {{'weights'}, {'array of weights, [npts,1]'}};
        {{'param_coords'}, {'parametric coordinates of quadrature points, array [npts,dim]'}}};
    return;
else
    prop_name=varargin{1};
    switch prop_name
        case {'npts','npts_per_gcell'}
            val = length(self.weights);
        case 'weights'
            val = self.weights;
        case 'param_coords'
            val = self.param_coords;
        otherwise
            error(['Unknown property name ''' prop_name '''!']);
    end
end
return;

