function magnitudefield = magnitude(self,varargin)
% Compute the magnitudes of the field values.
%
% function magnitudefield = magnitudes(self,varargin)
%
%   The new field on output holds for each node the magnitude
%   of the nodal degree of freedom in the appropriate norm (default is 2-norm).
%
    which_norm = 2;
    if (nargin > 2)
        which_norm = varargin{1};
    end
    wv = self.values;
    wmag =zeros(size(wv, 1),1);
    for j=1:size(wv, 1)
        wmag(j) =  norm(wv(j,:),which_norm);
    end
    magnitudefield = nodal_field(struct ('name',['mag'], 'dim', 1, 'data',wmag));
end


