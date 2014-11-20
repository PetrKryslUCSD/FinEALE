function [h8fens,h8fes] = H8_extrude_Q4(fens,fes,nLayers,extrusionh)
% Extrude a mesh of quadrilaterals into a mesh of hexahedra (H8).
%
% function [fens,fes] = H8_extrude_Q4(fens,fes,nLayers,extrusionh)
%
% Inputs:
% fens,fes = nodes and cells of the quadrilateral mesh,
% nLayers= number of layers,
% extrusionh= handle to function that computes the extrusion vector
%   nxyz=extrusionh(xyz,k)
% 
% Output:
% fens= node set which incorporates the nodes on input plus the nodes of
% the extruded elements. Note that no merging (fusing) of nodes is performed.
% h8fes = H8 extruded geometric cells
%
% Examples: 
%     L= 0.3; % in-plane dimension
%     W = 0.3; % in-plane dimension
%     a= 0.15; % hole radius
%     H = 0.03; % thickness of the plate
%     nL=5;nH=1;nW=5;na=7;
%     [fens,fes]=Q4_elliphole(a,a,L,W,nL,nW,na,[]);
%     [fens,fes] = H8_extrude_Q4(fens,fes,nH,@(x,i)([x,0]+[0,0,H*i]));
%     drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
%
    id = (1:count(fens))';
    cn=connected_nodes(fes);
    id(cn)=(1:length(cn))';
    q4fes =update_conn(fes,id);
    xyz=fens.xyz;
    q4fens =fenode_set(struct('xyz',xyz(cn,:)));
    
    [h8fens,h8fes]= do_extrude(q4fens,q4fes,nLayers,extrusionh);
end

function [efens,efes]= do_extrude(fens,fes,nLayers,extrusionh)
    cn=connected_nodes(fes);
    xyz1=fens.xyz;
    nn1 =count(fens);
    nnt=nn1*nLayers;
    ngc=count(fes)*nLayers;
    hconn=zeros(ngc,8);
    xyz1=fens.xyz;
    xyz =zeros(size(xyz1,1)*(nLayers+1),3);
    for j=1:length(cn)
        xyz(cn(j),:) =extrusionh(xyz1(cn(j),:),0);
    end
    for k=1:nLayers
        for j=1:length(cn)
            f=cn(j)+k*nn1;
            xyz(f,:) =extrusionh(xyz1(cn(j),:),k);
        end
    end
    
    conn = fes.conn;
    gc=1;
    for k=1:nLayers
        for i=1:count(fes)
            hconn(gc,:)=[conn(i,:)+(k-1)*nn1, conn(i,:)+k*nn1];
            gc=gc+1;
        end
    end
    efes =fe_set_H8(struct('conn',hconn));
    efens =fenode_set(struct('xyz',xyz));
end