% Standard NASTRAN colormap.
    %
    %     b = nastrancolors(This_is_ignored)
    %
    %   cadcolors returns an M-by-3 matrix containing a "nastrancolors" colormap.
    %   If no figure exists, MATLAB creates one.
    %
    %   For example, to reset the colormap of the current figure:
    %
    %             colormap(nastrancolors)
    %
    %   See also HSV, GRAY, HOT, COOL, COPPER, PINK, FLAG,
    %   bwr, cadcolors, cadcolors2, nastrancolors,
    %   COLORMAP, RGBPLOT.
    
function m=nastrancolors(This_is_ignored)
m = [255,255,255;
    0,239,255;
    0,159,255;
    0,64,255;
    0,0,255;
    0,255,0;
    0,170,0;
    8,85,9;
    255,170,255;
    255,85,255;
    255,35,255;
    255,255,0;
    255,170,0;
    255,85,0;
    255,0,0;
    ]/255;
