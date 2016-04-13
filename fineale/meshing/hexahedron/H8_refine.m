function [fens,fes] = H8_refine(fens,fes)
% Refine a mesh of H8 hexahedrals by octasection.
%
% function [fens,fes] = H8_refine(fens,fes)
%
% Arguments and
% Output:
% fens= finite element node set
% fes = finite element set
%
% Examples: 
%
%     xyz = [3, 1, 6; -5, 2, 1];
%     [fens,fes] = H8_hexahedron(xyz,1,2,3);
%     drawmesh({fens,fes},'fes','facecolor','red');
%     [fens,fes] = H8_refine(fens,fes);
%     figure;
%     drawmesh({fens,fes},'fes','facecolor','m');

    [fens,fes] = H8_to_H27(fens,fes);
    conn=fes.conn; label=fes.label;
    nconn=zeros(8*size(conn,1),8);
    nlabel=[];
    if (~isempty(label))
        nlabel=reshape(repmat(label',8,1),[],1);
    end
    nc=1;
    for i= 1:size(conn,1)
        conn27=conn(i,:);
        nconn(nc,:) =conn27([1,9,21,12,17,22,27,25]);        nc= nc+ 1;
        nconn(nc,:) =conn27([9,2,10,21,22,18,23,27]);        nc= nc+ 1;
        nconn(nc,:) =conn27([21,10,3,11,27,23,19,24]);        nc= nc+ 1;
        nconn(nc,:) =conn27([12,21,11,4,25,27,24,20]);        nc= nc+ 1;
        nconn(nc,:) =conn27([17,22,27,25,5,13,26,16]);        nc= nc+ 1;
        nconn(nc,:) =conn27([22,18,23,27,13,6,14,26]);        nc= nc+ 1;
        nconn(nc,:) =conn27([27,23,19,24,26,14,7,15]);        nc= nc+ 1;
        nconn(nc,:) =conn27([25,27,24,20,16,26,15,8]);        nc= nc+ 1;
    end
    fes=fe_set_H8(struct('conn',nconn));
    if (~isempty(label))
        fes.label=nlabel;
    end
end
