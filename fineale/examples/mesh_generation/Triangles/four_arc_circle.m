% Mesh of a circle composed of four arcs.

h=2.8;
[fens,fes,groups,edge_fes,edge_groups]=targe2_mesher({...
    'curve 1 arc -12 0 0 -12 Center 0 0 ',...
    'curve 2 arc 0 -12 12 0 Center 0 0 ',...
    'curve 3 arc 12 0 0 12 Center 0 0 ',...
    'curve 4 arc 0 12 -12 0 Center 0 0 ',...
    ['subregion 1  property 1 boundary 1 2 3 4'],...
    ['m-ctl-point constant ' num2str(h)]
    }, 1.0);
drawmesh({fens,fes},'fes','shrink', 0.8);
view(2)