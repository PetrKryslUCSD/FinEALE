% Standard CAD 12-color colormap.
    %
    %     b = cadcolors2(This_is_ignored)
    %
    %   cadcolors2 returns an 12-by-3 matrix containing a "cadcolors2" colormap.
    %   If no figure exists, MATLAB creates one.
    %
    %   For example, to reset the colormap of the current figure:
    %
    %             colormap(cadcolors2)
    %
    %   See also HSV, GRAY, HOT, COOL, COPPER, PINK, FLAG,
    %   bwr, cadcolors, cadcolors2, nastrancolors,
    %   COLORMAP, RGBPLOT.
    
function m=cadcolors2(This_is_ignored)
m = [0         0    1.0000
         0    0.6667    1.0000
         0    1.0000    1.0000
    0.3333    1.0000    0.6667
    0.5    1.0000    0.6667
    0.6667    1.0000   0.5 
    0.6667    1.0000    0.3333
    1.0000    1.0000         0
    1.0000    0.6667         0
    1.0000    0.        0
    1.0000    0.        0];
