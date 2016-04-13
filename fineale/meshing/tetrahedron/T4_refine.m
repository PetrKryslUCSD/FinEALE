function [fens,fes] = T4_refine(fens,fes)
% Refine a mesh of 4-node tetrahedra by octasection.
%
% function [fens,fes] = T4_refine(fens,fes)
%
% Input/Output:
% fens= finite element node set
% fes = finite element set
%
% Examples:
%
% [fens,fes] = T4_cylinderdel(3.1,1.7,2,2);
% figure; drawmesh({fens,fes},'fes','facecolor','m'); hold on
% [fens,fes] = T4_refine(fens,fes);
% figure; drawmesh({fens,fes},'fes','facecolor','b'); hold on

    [fens,fes] = T4_to_T10(fens,fes);
    conn=fes.conn; label=fes.label;
    nconn=zeros(8*size(conn,1),4);
    nlabel=[];
    if (~isempty(label))
        nlabel=reshape(repmat(label',8,1),[],1);
    end
    
    nc=1;
    for i= 1:size(conn,1)
        conn10=conn(i,:);
        nconn(nc,:) =conn10([1,5,7,8]);        nc= nc+ 1;
        nconn(nc,:) =conn10([2,6,5,9]);        nc= nc+ 1;
        nconn(nc,:) =conn10([3,7,6,10]);        nc= nc+ 1;
        nconn(nc,:) =conn10([8,9,10,4]);        nc= nc+ 1;
        nconn(nc,:) =conn10([6,7,5,9]);        nc= nc+ 1;
        nconn(nc,:) =conn10([7,6,10,9]);        nc= nc+ 1;
        nconn(nc,:) =conn10([9,7,5,8]);        nc= nc+ 1;
        nconn(nc,:) =conn10([7,9,10,8]);        nc= nc+ 1;
    end
    fes=fe_set_T4(struct('conn',nconn));
    if (~isempty(label))
        fes.label=nlabel;
    end
    
end
