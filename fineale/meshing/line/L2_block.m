function [fens,fes] = L2_block(Length,nL,options)
% Mesh of a line segment located in a given range, L2 Elements
%
% function [fens,fes] = L2_block(Length,nL,options)
%
% Interval x=<0,Length>
% Divided into nL elements of the L2 type.
%
% Examples: 
%     [fens,fes] = L2_block(5.0,4,struct('other_dimension',.25));
%     drawmesh({fens,fes},'nodes','fes','facecolor','none', 'linewidth',3); hold on
%
%
% See also: L2_blockx, fe_set_L2
%
    [fens,fes] = L2_blockx(linspace(0,Length,nL+1)', options);
end
