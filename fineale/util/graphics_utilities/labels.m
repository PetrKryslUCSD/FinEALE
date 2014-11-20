function labels (lx,ly,lz)
% Set labels on one, two or all three axes.
%
% function labels (lx,ly,lz)
%
% If lx is supplied as empty,
%
% labels ([])
%
% or if the function is called without arguments
%
% the default labels (X, Y, Z) are used in all three directions.
%
% Note the interpreter is latex by default.
    if (~exist ('lx','var'))
        lx='X'; ly='Y'; lz='Z';
    end
    if (~exist ('ly','var'))
        ly = [];
    end
    if (~exist ('lz','var'))
        lz = [];
    end
    if (~isempty(lx))
        l=xlabel(lx);
        set(l,'interpreter','latex')
    end
    if (~isempty(ly))
        l=ylabel(ly);
        set(l,'interpreter','latex')
    end
    if (~isempty(lz))
        l=zlabel(lz);
        set(l,'interpreter','latex')
    end
end
function s1 = strim(s)
%STRIM(S) strip the trailing and leading blanks of a string.
%         STRIM(S)is an extension of the MATLAB function DEBLANK(S)
%         and removes both leading and trailing blanks from the string S.
%
% Date 11-08-98
% Written by Stefan Baunack. Please send any bug
% reports or other comments to: s.baunack@ifw-dresden.de.

if ~isstr(s)
	error('Input must be a string.')
end
s1 = deblank(s);
s1 = fliplr(s1);
s1 = deblank(s1);
s1 = fliplr(s1);
end