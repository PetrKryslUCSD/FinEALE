% Get property from the specified object and return the value.
%   Call as:
%      get(obj), or
%      get(obj,prop_name)
%   where
%      obj       - object
%      prop_name - property name
%
%   The first form returns the list of property names as a string,
%   the second form returns the value of the named property.
%
function val = get(self,varargin)
    if (nargin == 1)
        val = {{{'sigma_y'}, {'yield stress, scalar'}}
            };
        val = cat(1, val, get(self.property_linel_iso));
        return;
    else
        prop_name=varargin{1};
        switch prop_name
            case 'sigma_y'
                val = self.sigma_y;
            case 'Hi'
                val = self.Hi;
            case 'Hk'
                val = self.Hk;
            case 'Hn'
                val = self.Hn;
            otherwise
                val = get(self.property_linel_iso, varargin{:});
        end
    end
end
