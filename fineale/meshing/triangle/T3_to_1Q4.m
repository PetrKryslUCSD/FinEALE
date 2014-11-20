function [fens,fes] = T3_to_1Q4(fens,fes,options)
% Convert a mesh of triangles T3 to one quadrilateral Q4 each.
%
% function [fens,fes] = T3_to_1Q4(fens,fes,options)
%
% options =attributes recognized by the constructor fe_set_Q4
    if ~isstruct(options)
        options = struct('conn',[]);
    end
    % construct new geometry cells
    nconns=zeros(count(fes),4);
    conns = fes.conn;
    nc =zeros(length(fens),1);
    for i= 1:count(fes)
        conn=conns(i,:);
        nc(conn) =nc(conn)+1;
    end
    gc=1;
    for i= 1:size(conns,1)
        conn=conns(i,:);
        nc1 =nc(conn);
        [C,I] = max(nc1); 
        if (C<2)
            error(' Cannot avoid having only collapsed edges at node');
        end
        nconns(gc,:) =conn([1:I,I:3]);
        gc=gc+1;
        nc(conn(I)) =nc(conn(I))-1;
    end
    options.conn=nconns;
    fes=fe_set_Q4(options);
end
