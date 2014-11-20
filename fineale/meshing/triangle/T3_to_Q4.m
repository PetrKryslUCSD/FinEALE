function [fens,fes] = T3_to_Q4(fens,fes,options)
% Convert a mesh of triangles T3 to three quadrilateral Q4 each.
%
% function [fens,fes] = T3_to_Q4(fens,fes,options)
%
% options =attributes recognized by the constructor fe_set_Q4
    [fens,fes] = T3_to_T6(fens,fes,options);
    if ~isstruct(options)
        options = struct('conn',[]);
    end
    onfens=count(fens);
    xyz1=fens.xyz;
    id1=(1:count(fens))';
    id=zeros(count(fens)+count(fes),1);
    id(1:count(fens),1)=id1;
    xyz=zeros(count(fens)+count(fes),3);
    xyz(1:count(fens),1:size(xyz1,2))=xyz1;
    % construct new geometry cells
    nconns=zeros(3*count(fes),4);
    conns = fes.conn;
    nc=onfens+1;
    gc=1;
    for i= 1:count(fes)
        conn=conns(i,:);
        id(nc)=nc;
        xyz(nc,:)=mean(xyz(conn(1:3),:));
        nconns(gc,:) =[conn(1) conn(4) nc conn(6)];
        gc=gc+1;
        nconns(gc,:) =[conn(2) conn(5) nc conn(4)];
        gc=gc+1;
        nconns(gc,:) =[conn(3) conn(6) nc conn(5)];
        gc=gc+1;
        nc= nc+ 1;
    end
    options.conn=nconns;
    fes=fe_set_Q4(options);
    fens.xyz=xyz;
end
