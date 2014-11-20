function  [fens, fes] = L2_blockx(xyz, options)
% Graded mesh of a 1-D block, L2 finite elements. 
%
% Syntax
%
% function [fens,fes] = L2_blockx(xyz,options)
%
% Input
%
% xyz - (n x p) matrix containing the locations of the individual nodes.
%
% options - structure with fields recognized by the constructor of the
%   L2 finite element
%
% Example
% 
%     [fens,fes] = L2_blockx((2:1:7)'.^3, struct('other_dimension', 0.1));
%     drawmesh({fens,fes},'nodes','fes','facecolor','none', 'linewidth',3); hold on
%
%     a=pi/180*(0:1:6)'/6*135;
%     [fens,fes] = L2_blockx(5.1*[sin(a),cos(a)], struct('other_dimension', 0.1));
%     drawmesh({fens,fes},'nodes','fes','facecolor','none', 'linewidth',3); hold on
% See also: L2_block, fe_set_L2

    if ~isstruct(options)
        other_dimension = options; 
        clear options;
        options.other_dimension = other_dimension;
    end
    
    ncells = length(xyz) - 1;
    nnodes = (ncells+1);

    % create the nodes
    fens = fenode_set(struct('xyz', xyz));
    
    % create the cells
    options.conn = [(1:ncells)', (2:ncells+1)'];
    fes = fe_set_L2(options);
    
end