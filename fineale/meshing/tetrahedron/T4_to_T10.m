function [fens,fes] = T4_to_T10(fens,fes)
% Convert a mesh of Tetrahedron T4 (four-node) to Tetrahedron T10.
%
% function [fens,fes] = T4_to_T10(fens,fes)
%
% Examples: 
% [fens,fes] = T4_sphere(3.1,1);
% [fens,fes] = T4_to_T10(fens,fes);
% fens= onto_sphere(fens,3.1,connected_nodes(mesh_boundary(fes,[])));
% figure; drawmesh({fens,fes},'fes','facecolor','y'); hold on
%
conn =fes.conn;
    label =fes.label;
    nedges=6;
    ec = [1, 2; 2, 3; 3, 1; 4, 1; 4, 2; 4, 3];
    % make a search structure for edges
    edges={};
    for i= 1:size(conn,1)
        for J = 1:nedges
            ev=conn(i,ec(J,:));
            anchor=min(ev);
            if length(edges)<anchor
                edges{anchor}=[];
            end
            edges{anchor}=unique([edges{anchor} max(ev)]);
        end
    end
    % now generate new node number for each edge
    xyz =fens.xyz;
    n=size(xyz,1);
    nodes=edges;
    for i= 1:length(edges)
        e=edges{i};
        for J = 1:length(e)
            n=n+1;
            e(J)=n;
        end
        nodes{i}=e;
    end
    xyz(size(xyz,1)+1:n,1:3)=0;
    % calculate the locations of the new nodes
    % and construct the new nodes
    for i= 1:length(edges)
        e=edges{i};
        n=nodes{i};
        for J = 1:length(e)
            id(n(J))=n(J);
            xyz(n(J),:)=(xyz(i,:)+xyz(e(J),:))/2;
        end
    end
    nfens=length(id);% number of nodes in the original mesh plus number of the edge nodes
    fens=fenode_set(struct('xyz',xyz));
    % construct new geometry cells
    nconn=zeros(size(conn,1),10);
    nc=1;
    for i= 1:size(conn,1)
        econn=zeros(1,nedges);
        for J = 1:nedges
            ev=conn(i,ec(J,:));
            anchor=min(ev);
            e=edges{anchor};
            n=nodes{anchor};
            econn(J)=n(find(e==max(ev)));
        end
        nconn(nc,:) =[conn(i,:) econn];
        nc= nc+ 1;
    end
    fes=fe_set_T10(struct('conn',nconn,'label',label));
    return;
end
