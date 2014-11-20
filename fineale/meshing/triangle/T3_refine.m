function [fens,fes] = T3_refine(fens,fes,options)
% Refine a mesh of 3-node tetrahedra by quadrisection.
%
% function [fens,fes] = T3_refine(fens,fes,options)
%
% Input/Output:
% fens= finite element node set
% fes = finite element set
% options =  struct recognized by the constructor of T3
%
% Examples:
%
% [fens,fes] = T4_cylinderdel(3.1,1.7,2,2);
% figure; drawmesh({fens,fes},'fes','facecolor','m'); hold on
% [fens,fes] = T4_refine(fens,fes);
% figure; drawmesh({fens,fes},'fes','facecolor','b'); hold on

    [fens,fes] = T3_to_T6(fens,fes,options);
    conn=fes.conn;
    nconn=zeros(4*size(conn,1),3);
    nc=1;
    for i= 1:size(conn,1)
        c=conn(i,:);
        nconn(nc,:) =c([1,4,6]);        nc= nc+ 1;
        nconn(nc,:) =c([2,5,4]);        nc= nc+ 1;
        nconn(nc,:) =c([3,6,5]);        nc= nc+ 1;
        nconn(nc,:) =c([4,5,6]);        nc= nc+ 1;
    end
    fes=fe_set_T3(struct('conn',nconn));
end
