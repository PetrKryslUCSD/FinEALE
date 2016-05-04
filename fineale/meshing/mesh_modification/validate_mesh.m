function validate_mesh(fens, fes)
% Validate finite element mesh.
%
% function validate_mesh(fens, fes)
%
% fes= either a finite element set, or a cell array of finite element sets
%
% A finite element mesh given by the node set and the finite element set
% is validated by checking the sanity of the numbering:
% -- the node numbers need to be positive and in serial order
% -- the fe connectivity needs to refer to valid nodes
% -- the finite element nodes need to be connected to at least one finite
% element
%
% An error is reported as soon as it is detected.
% 
    totnfens=count(fens);
    
    %     For uniformity of treatment we will convert the finite elements
    %     into a cell array of finite element sets
    if (~iscell(fes))
        fes={fes};
    end
    
    connected=[];
    for   fi=1:length(fes)
        conn=fes{fi}.conn;
        
        %     sconn =sort (conn,2);
        %     uconn =unique(sconn,'rows');
        %     if (size(conn,1)~=size(uconn,1))
        %         error (['Repeated connectivity: ' num2str(size(uconn,1)) ' of ' num2str(size(conn,1)) ' unique'])
        %     end
        
        for i=1:size(conn,1)
            if (max(conn(i,:))>totnfens)
                error (['Wrong connectivity (refers to nonexistent node) ' num2str(conn(i,:))])
            end
            if (min(conn(i,:))<1)
                error (['Wrong connectivity (refers to nonexistent node) ' num2str(conn(i,:))])
            end
            if (length(unique(conn(i,:)))~=length(conn(i,:)))
                error (['Wrong connectivity (multiply referenced node) ' num2str(conn(i,:))])
            end
        end
        
        connected1 = find_unconn_fens(fens, fes{fi});
        if (isempty( connected ))
            connected= connected1;
        end
        connected= connected+connected1;
       
    end% for all sets
    connected(connected>0)=1;
    
    %     Check whether any nodes are without any connections by elements: if so they need to be removed.
    if (any(connected==0))
        error (['Unconnected nodes present: ' num2str(sum(connected==0)) ' total'])
    end
    
end
