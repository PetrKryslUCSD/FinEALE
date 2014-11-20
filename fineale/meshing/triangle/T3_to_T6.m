function [fens,fes] = T3_to_T6(fens,fes,options)
% Convert a mesh of triangle T3 (three-node) to triangle T6.
%
% function [fens,fes] = T3_to_T6(fens,fes,options)
%
% options =attributes recognized by the constructor fe_set_T6
    if ~isstruct(options)
        options = struct('conn',[]);
    end
    nedges=3;
    ec = [1, 2; 2, 3; 3, 1];
    % make a search structure for edges
    edges={};
    conns = fes.conn;
    label =fes.label;
    for i= 1:count(fes)
        conn=conns(i,:);
        for J = 1:nedges
            ev=conn(ec(J,:));
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
    xyz=fens.xyz;
    xyz =zeros(n,size(xyz,2));
    xyz(1:count(fens),:) = fens.xyz;
    id =zeros(n,1); 
    id = (1:count(fens))';
    % calculate the locations of the new nodes
    % and construct the new nodes
    for i= 1:length(edges)
        e=edges{i};
        n=nodes{i};
        for J = 1:length(e)
            xyz(n(J),:)=(xyz(i,:)+xyz(e(J),:))/2;
            id(n(J),:)=n(J);
        end
    end
    nfens=count(fens);% number of nodes in the original mesh plus number of the edge nodes
    fens.xyz=xyz;
    % construct new geometry cells
    nconns=zeros(count(fes),6);
    nc=1;
    for i= 1:count(fes)
        conn = conns(i,:);
        econn=zeros(1,nedges);
        for J = 1:nedges
            ev=conn(ec(J,:));
            anchor=min(ev);
            e=edges{anchor};
            n=nodes{anchor};
            econn(J)=n(find(e==max(ev)));
        end
        nconns(nc,:) =[conn econn];
        nc= nc+ 1;
    end
    Parameters =options;
    Parameters.conn=nconns;
    Parameters.label=label;
    fes=fe_set_T6(Parameters) ;
end
