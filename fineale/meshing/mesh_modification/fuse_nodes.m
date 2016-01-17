function [fens,new_indexes_of_fens1_nodes] = fuse_nodes(fens1, fens2, tolerance)
    % Fuse together nodes from to node sets.
    %
    % function [fens,new_indexes_of_fens1_nodes] = fuse_nodes(fens1, fens2, tolerance)
    %
    % Fuse two node sets. If necessary, by gluing together nodes located within tolerance of each other.
    % The two node sets, fens1 and fens2,  are fused together by
    % merging the nodes that fall within a box of size "tolerance".
    % The merged node set, fens, and the new  indexes of the nodes
    % in the set fens1 are returned.
    %
    % The set fens2 will be included unchanged, in the same order,
    % in the node set fens.
    %
    % The indexes of the node set fens1 will have changed.
    %
    % Example: 
    % After the call to this function we have
    % k=new_indexes_of_fens1_nodes(j) is the node in the node set fens which
    % used to be node j in node set fens1.
    % The finite element set connectivity that used to refer to fens1
    % needs to be updated to refer to the same nodes in  the set fens as
    %     fes = update_conn(fes ,new_indexes_of_fens1_nodes);
    %
    % See also: merge_nodes, update_conn
    %
    
    
    xyz1 = fens1.xyz;
    id1 =(1:size(xyz1,1))';
    dim =size(xyz1,2);
    xyz2 = fens2.xyz;
    id2 =(1:size(xyz2,1))';
    c1=ones(size(xyz2,1),1);
    % Mark nodes from the first array that are duplicated in the second
    if (tolerance>0)% should we attempt to merge nodes?
        Box2= inflate_box(bounding_box(xyz2), tolerance);
        for i=1:count(fens1)
            XYZ =xyz1(i,:);
            % Note  that we are looking for  distances  of this node to nodes in the OTHER node set
            if (in_box (Box2,XYZ))% This check makes this run much faster
                xyzd=abs(xyz2-c1*XYZ);% find the distances along  coordinate directions
                jx=find(sum(xyzd')'<tolerance);
                if (~isempty(jx))
                    minn=min(jx);
                    id1(i) =-minn;
                end
            end
        end
        %                 for i=1:count(fens1)
        %                         XYZ =xyz1(i,:);
        %                         % Note  that we are looking for  distances  of this node to nodes in the OTHER node set
        %                         xyzd=abs(xyz2-c1*XYZ);% find the distances along  coordinate directions
        %                         jx=find(sum(xyzd')'<tolerance);
        %                         if (~isempty(jx))
        %                             minn=min(jx);
        %                             id1(i) =-minn;
        %                         end
        %                     end
    end
    % Generate  fused arrays of the nodes
    xyzm = zeros(count(fens1)+count(fens2),dim);
    xyzm(1:count(fens2),1:size(xyz2,2))=xyz2;% fens2 are there without change
    idm = zeros(count(fens1)+count(fens2),1);
    idm(1:count(fens2))=(1:count(fens2));
    mid=count(fens2)+1;
    for i=1:count(fens1) % and then we pick only non-duplicated fens1
        if id1(i)>0
            id1(i)=mid;
            idm(mid)=mid;
            xyzm(mid,:)=xyz1(i,:);
            mid=mid+1;
        else
            id1(i)=id2(-id1(i));
        end
    end
    nfens =mid-1;
    xyzm =xyzm(1:nfens,:);
    idm =idm(1:nfens,:);
    
    % Create the fused Node set
    fens =fenode_set(struct('xyz',xyzm));
    % The Node set 1 numbering will change
    new_indexes_of_fens1_nodes=id1;
    % The node set 2 numbering stays the same
end

