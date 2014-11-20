% Test generation of triangular mesh on a circle, visualization test.

% The circle circumference consists of two arcs, which are given centers
% slightly offset from the true center of the circle.  The mesh generator
% needs this to make sense of the orientation of the arcs.
h=3.2;
[fens,fes,groups,edge_fes,edge_groups]=targe2_mesher({...
    'curve 1 arc -12 0 12 0 Center 0 0.000001 ',...
    'curve 2 arc 12 0 -12 0 Center 0 -0.000001',...
    ['subregion 1  property 1 boundary  1  2 '],...
    ['m-ctl-point constant ' num2str(h)]
    }, 1.0);
drawmesh({fens,fes},'fes','shrink', 0.8);
view(2);pause(1);

drawmesh({fens,fes},'fes','shrink',1);
view(2);pause(1);

drawmesh({fens,fes},'fes','facecolor','r');
view(2);pause(1);

drawmesh({fens,fes},'fes','facecolor','r','edgecolor','none');
view(2);pause(1);

drawmesh({fens,fes},'fes','facecolor','none','edgecolor','none');
view(2);pause(1);

drawmesh({fens,fes},'fes','facecolor','y','shrink', 0.8,'edgecolor','none');
view(2);pause(1);
