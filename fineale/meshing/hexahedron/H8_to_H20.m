function [fens,fes] = H8_to_H20(fens,fes)
% Convert a mesh of hexahedra H8 to hexahedra H20.
%
% function [fens,fes] = H8_to_H20(fens,fes)
%
% Arguments and
% Output:
% fens= finite element node set
% fes = finite element set
%
    nedges=12;
    ec = [1, 2; 2, 3; 3, 4; 4, 1; 5, 6; 6, 7; 7, 8; 8, 5; 1, 5; 2, 6; 3, 7; 4, 8;];
    % make a search structure for edges
    edges={};
    conn = fes.conn;
    label = fes.label;
    for i= 1:count(fes)
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
    n=count(fens);
    nodes=edges;
    for i= 1:length(edges)
        e=edges{i};
        for J = 1:length(e)
            n=n+1;
            e(J)=n;
        end
        nodes{i}=e;
    end
    nfens=n;% number of nodes in the original mesh plus number of the edge nodes
    xyz =zeros(n,3);
    xyz(1:count(fens),:) = fens.xyz;
    % calculate the locations of the new nodes
    % and construct the new nodes
    for i= 1:length(edges)
        e=edges{i};
        n=nodes{i};
        for J = 1:length(e)
            xyz(n(J),:)=(xyz(i,:)+xyz(e(J),:))/2;
        end
    end
    fens.xyz=xyz;
     % construct new geometry cells
    nconns=zeros(count(fes),20);
    conns = fes.conn;
    for i= 1:count(fes)
        conn =conns(i,:);
        econn=zeros(1,nedges);
        for J = 1:nedges
            ev=conn(ec(J,:));
            anchor=min(ev);
            e=edges{anchor};
            n=nodes{anchor};
            econn(J)=n(e==max(ev));
        end
        nconns(i,:)=[conn econn];
    end
    fes=fe_set_H20(struct('conn',nconns,'label',label)) ;
    return;
end
