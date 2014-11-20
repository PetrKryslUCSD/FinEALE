% Standard ABAQUS 12-color colormap.
    %
    %     b = abaquscolors(This_is_ignored)
    %
    %   abaquscolors returns an 12-by-3 matrix containing a "abaquscolors" colormap.
    %   If no figure exists, MATLAB creates one.
    %
    %   For example, to reset the colormap of the current figure:
    %
    %             colormap(abaquscolors)
    %
    %   See also HSV, GRAY, HOT, COOL, COPPER, PINK, FLAG,
    %   bwr, cadcolors, cadcolors2, nastrancolors,
    %   COLORMAP, RGBPLOT.
    
function m=abaquscolors(This_is_ignored)
m = [0,0,0
    0         0    1.0000
         0    0.6667    1.0000
         0    1.0000    1.0000
    0.3333    1.0000    0.6667
    0.5    1.0000    0.6667
    0.6667    1.0000   0.5 
    0.6667    1.0000    0.3333
    1.0000    1.0000         0
    1.0000    0.6667         0
    1.0000    0.        0
    1.0000    0.        0
    00.75, 0.75, 0.75];
