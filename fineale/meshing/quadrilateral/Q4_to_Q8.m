function [fens,fes] = Q4_to_Q8(fens,fes,options)
% Convert a mesh of quadrilateral Q4 to quadrilateral Q8.
%
% function [fens,fes] = Q4_to_Q8(fens,fes,options)
%
% options =attributes recognized by the constructor fe_set_Q8
%
% Examples: 
%     R=8.5;
%     [fens,fes]=Q4_sphere(R,1,1.0);
%     [fens,fes] = Q4_to_Q8(fens,fes,[]);
%     fens= onto_sphere(fens,R,[]);
%     drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
%
% See also: fe_set_Q8

if (isempty(options))
	options=1.0;
end
if ~isstruct(options)
	other_dimension = options; clear options;
	options.other_dimension = other_dimension;
end
    nedges=4;
    ec = [1, 2; 2, 3; 3, 4; 4, 1];
     conns = fes.conn;
       % make a search structure for edges
    edges={};
    for i= 1:count(fes)
        conn = conns(i,:);
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
    xyz1=fens.xyz;
    xyz =zeros(n,size(xyz1,2));
    xyz(1:size(xyz1,1),:)=xyz1;
    % calculate the locations of the new nodes
    % and construct the new nodes
    for i= 1:length(edges)
        e=edges{i};
        n=nodes{i};
        for J = 1:length(e)
            xyz(n(J),:)=(xyz(i,:)+xyz(e(J),:))/2;
        end
    end
    nfens=size(xyz,2);% number of nodes in the original mesh plus number of the edge nodes
    % construct new geometry cells
    nconns =zeros(count(fes),8);
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
    Parameters=options;
    Parameters.other_dimension=fes.other_dimension;
    Parameters.conn=nconns;
    fes =fe_set_Q8(Parameters);
    fens =fenode_set(struct('xyz',xyz));
end
