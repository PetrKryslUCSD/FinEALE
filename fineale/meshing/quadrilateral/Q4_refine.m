% Refine a mesh of quadrilaterals by bisection
%
% function [fens,fes] = Q4_refine(fens,fes)
%
% Examples: 
% [fens,fes] = Q4_quadrilateral([-1,-1;2,-2;3,3;-1,1],2,3,[]);
% [fens,fes] = Q4_refine(fens,fes);
% drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on

function [fens,fes] = Q4_refine(fens,fes)
    nedges=4;
    nsubcells=4;
    % make a search structure for edges
    edges={};
    conn = fes.conn;
    for i= 1:size(conn,1)
        for J = 1:nedges
            ev=conn(i,[J, mod(J,nedges)+1]);
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
    xyz =fens.xyz;
    xyz(count(fens)+1:n,:)=zeros(n-count(fens),size(xyz,2));
    % calculate the locations of the new nodes
    % and construct the new nodes
    for i= 1:length(edges)
        e=edges{i};
        n=nodes{i};
        for J = 1:length(e)
            xyz(n(J),:)=mean(xyz([i,e(J)],:));
        end
    end
    nfens=size(xyz,1);% number of nodes in the original mesh plus number of the edge nodes
    nifens=count(fes);% number of internal nodes
    xyz(nfens+1:nfens+nifens,:)=zeros(nifens,size(xyz,2));
    nfens=size(xyz,1);
    % construct new geometry cells
    nconn =zeros(nsubcells*count(fes),4);
    nc=1;
    for i= 1:count(fes)
        econn=zeros(1,nedges);
        for J = 1:nedges
            ev=conn(i,[J, mod(J,nedges)+1]);
            anchor=min(ev);
            e=edges{anchor};
            n=nodes{anchor};
            econn(J)=n(find(e==max(ev)));
        end
        xyzm=mean(xyz(conn(i,:),:));
        inn=nfens-nifens+i;
        xyz(inn,:)=xyzm;
        nconn(nc,:) =[conn(i,1) econn(1) inn econn(4)];
        nc= nc+ 1;
        nconn(nc,:) =[conn(i,2) econn(2) inn econn(1)];
        nc= nc+ 1;
        nconn(nc,:) =[conn(i,3) econn(3) inn econn(2)];
        nc= nc+ 1;
        nconn(nc,:) =[conn(i,4) econn(4) inn econn(3)];
        nc= nc+ 1;
    end
    Thickness =fes.eval_other_dimension;
    fes =fe_set_Q4(struct('conn',nconn,...
        'other_dimension',Thickness,'axisymm', fes.axisymm));
    fens =fenode_set(struct('xyz',xyz));
end
