function [fens,fes] = Q4_L2x2
% L-shaped  domain mesh.
%
% function [fens,fes] = Q4_L2x2
%
% 2 x 2 quadrilateral element mesh of one quarter of a square domain with a square hole.
%
% Examples: 
% [fens,fes] = Q4_L2x2
% drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
    fensd = [f(1, [1/2   0]); ...
        f(2, [1   0]); ...
        f(3, [1   1/2]); ...
        f(4, [1/2   1/2]);...
        f(5, [1   1]);...
        f(6, [1/2   1]);...
        f(7, [0   1]);...
        f(8, [0   1/2]);...
        ];
    fes = [g(1, [1, 2, 3, 4], 1.0);
        g(1, [4, 3, 5, 6], 1.0);
        g(1, [4, 6, 7, 8], 1.0)];
    fens=fenode_set(struct('xyz',fensd(:,2:3)));
    options.conn =fes(:,2:5);
    fes =fe_set_Q4(options);
function v=f(ID,x)
    v=[ID,x];

function v=g(ID,conn, thickness)
    v=[ID,conn,thickness];
