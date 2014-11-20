function validate_mesh(fens, fes)
%  Validate finite element mesh.
%
% function validate_mesh(fens, fes)
%
% A finite element mesh given by the node set and the finite element set
% is validated by checking the sanity of the numbering:
% -- the node numbers need to be positive and in serial order
% -- the fe connectivity needs to refer to valid nodes
% 
% An error is reported as soon as it is detected.
% 
    nfens=count(fens);
    
    conn=fes.conn;
    
%     sconn =sort (conn,2);
%     uconn =unique(sconn,'rows');
%     if (size(conn,1)~=size(uconn,1))
%         error (['Repeated connectivity: ' num2str(size(uconn,1)) ' of ' num2str(size(conn,1)) ' unique'])
%     end
    
    for i=1:size(conn,1)
        if (max(conn(i,:))>nfens)
            error (['Wrong connectivity (refers to nonexistent node) ' num2str(conn(i,:))])
        end
        if (min(conn(i,:))<1)
            error (['Wrong connectivity (refers to nonexistent node) ' num2str(conn(i,:))])
        end
        if (length(unique(conn(i,:)))~=length(conn(i,:)))
            error (['Wrong connectivity (multiply referenced node) ' num2str(conn(i,:))])
        end
    end
end
