% Convert a mesh of quadrilateral Q4s into T3s using short diagonal.
%
% function [fens,fes] = Q4_to_T3_sd(fens,fes,options)
%
% The connectivity associated with a shorter diagonal in the quadrilateral
% is chosen for the split.
% 
% options =attributes recognized by the constructor fe_set_T3, and
% orientation = 'default' or 'alternate' chooses which diagonal is taken
%      for splitting
% 
% Example: 
%     [fens,fes] = Q4_quadrilateral([-1,-1;2,-2;3,3;-1,1],2,3,[]);
%     [fens,fes] = Q4_to_T3_sd(fens,fes,[]);
%     drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
%
% See also: Q4_to_T3
% 
function [fens,fes] = Q4_to_T3_sd(fens,fes,options)
    if ~isstruct(options)
        options = struct('id',1);
    end
    x=fens.xyz;
    connl1=(1:3);
    connl2=[1, 3, 4];
    connl1a=[1, 2, 4];
    connl2a=[3, 4, 2];
    nedges=4;
    nconns=zeros(2*count(fes),3);
    conns = fes.conn;
    nc=1;
    for i= 1:count(fes)
        conn = conns(i,:);
        c1=conn(connl1);
        c2=conn(connl2);
        el= [elens(x(c1,:)),elens(x(c2,:))];
        c1a=conn(connl1a);
        c2a=conn(connl2a);
        ela= [elens(x(c1a,:)),elens(x(c2a,:))];
        if (max(el)<max(ela))
            nconns(nc,:) =c1;
            nc= nc+ 1;
            nconns(nc,:) =c2;
            nc= nc+ 1;
        else
            nconns(nc,:) =c1a;
            nc= nc+ 1;
            nconns(nc,:) =c2a;
            nc= nc+ 1;
        end
    end
    options.conn=nconns;
    fes=fe_set_T3(options);
end

function el = elens(x)
    el=[norm(diff(x([1,2],:))),norm(diff(x([3,2],:))),norm(diff(x([1,3],:)))];
end
