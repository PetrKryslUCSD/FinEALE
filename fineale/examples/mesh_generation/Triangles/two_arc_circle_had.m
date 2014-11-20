% Graded mesh on a circle.
% The circle circumference consists of two arcs, which are given centers
% slightly offset from the true center of the circle.  The mesh generator
% needs this to make sense of the orientation of the arcs.
h=1.2;
[fens,fes,groups,edge_fes,edge_groups]=targe2_mesher({'regext -12.1 12.1 -12.1 12.1',...
    'curve 1 arc -12 0 12 0 Center 0 0.000001 ',...
    'curve 2 arc 12 0 -12 0 Center 0 -0.000001',...
    ['subregion 1  property 1 boundary  1  2 '],...
    ['m-ctl-point constant ' num2str(h)]
    }, 1.0);
drawmesh({fens,fes},'fes','shrink', 0.8);
view(2)

Mesh_size=@(x) 0.1+4.0*sin(0.5*norm(x))^2;

for k=1:3
    x=fens.xyz;
    h=zeros(count(fens),1);
    for j=1:count(fens)
        h(j)=Mesh_size(x(j,:));
    end
    conn=fes.conn;
    
    [fens,fes,groups,edge_fes,edge_groups] = targe2_mesher_adapt({'regext -12.1 12.1 -12.1 12.1',...
        'curve 1 arc -12 0 12 0 Center 0 0.000001 ',...
        'curve 2 arc 12 0 -12 0 Center 0 -0.000001',...
        ['subregion 1  property 1 boundary  1  2 '],...
        },conn,x,h, 1.0);
    drawmesh({fens,fes},'fes');
    view(2)
end
