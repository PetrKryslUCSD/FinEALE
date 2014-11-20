function m=cadcolors(This_is_ignored)
% Standard CAD 10-color colormap.
    %
    %     b = cadcolors(This_is_ignored)
    %
    %   cadcolors returns an M-by-3 matrix containing a "cadcolors" colormap.
    %   If no figure exists, MATLAB creates one.
    %
    %   For example, to reset the colormap of the current figure:
    %
    %             colormap(cadcolors)
    %
    %   See also HSV, GRAY, HOT, COOL, COPPER, PINK, FLAG,
    %   bwr, cadcolors, cadcolors2, nastrancolors,
    %   COLORMAP, RGBPLOT.
    
m = [ 0         0    1.0000
         0    0.6667    1.0000
         0    1.0000    1.0000
    0.3333    1.0000    0.6667
    0.5    1.0000    0.6667
    0.6667    1.0000   0.5 
    0.6667    1.0000    0.3333
    1.0000    1.0000         0
    1.0000    0.6667         0
    1.0000    0.        0];
