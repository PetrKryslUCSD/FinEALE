function [q4fens,q4fes] = Q4_extrude_L2(fens,fes,nLayers,extrusionh)
% Extrude a mesh of L2 line segments into a mesh of Q4 quadrilaterals.
%
% function [fens,fes] = Q4_extrude_L2(fens,fes,nLayers,extrusionh)
%
%
% Inputs:
% fens,fes = nodes and elements of the quadrilateral mesh,
% nLayers= number of layers,
% extrusionh= handle to function that computes the extrusion vector
%   nxyz=extrusionh(xyz,k)
% 
% Output:
% fens= node set which incorporates the nodes on input plus the nodes of
% the extruded elements. Note that no merging (fusing) of nodes is performed.
% q4fes = Q4 extruded finite elements
%
% Examples: 
%     [fens,fes] = L2_block(4.0,5, struct('other_dimension', 3.1));
%     [fens,fes] = Q4_extrude_L2(fens,fes,2,@(xyz, layer)[xyz(1),(xyz(1)+1/5)*layer]);
%     drawmesh({fens,fes},'fes','facecolor','red')
%
% See also: fe_set_Q4

    id = (1:count(fens))';
    cn=connected_nodes(fes);
    id(cn)=(1:length(cn))';
    L2fes =update_conn(fes,id);
    xyz=fens.xyz;
    L2fens =fenode_set(struct('xyz',xyz(cn,:)));
    
    [q4fens,q4fes]= do_extrude(L2fens,L2fes,nLayers,extrusionh);
end

function [efens,efes]= do_extrude(fens,fes,nLayers,extrusionh)
    nn1 =count(fens);
    nnt=nn1*nLayers;
    ngc=count(fes)*nLayers;
    xyz1=fens.xyz;
    x1=extrusionh(xyz1(1,:),0);
    xyz =zeros(size(xyz1,1)*(nLayers+1),size(x1,2));
    qconn=zeros(ngc,4);
    for j=1:nn1
        xyz(j,:) =extrusionh(xyz1(j,:),0);
    end
    for k=1:nLayers
        for j=1:nn1
            f=j+k*nn1;
            xyz(f,:) =extrusionh(xyz1(j,:),k);
        end
    end
    
    conn = fes.conn;
    gc=1;
    for k=1:nLayers
        for i=1:count(fes)
            qconn(gc,:)=[conn(i,:)+(k-1)*nn1, conn(i,end:-1:1)+k*nn1];
            gc=gc+1;
        end
    end
    efes =fe_set_Q4(struct('conn',qconn));
    efens =fenode_set(struct('xyz',xyz));
end
