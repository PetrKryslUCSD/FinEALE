function [fens,fes] = merge_nodes(fens, fes, tolerance)
    % Merge together  nodes of a single node set.
    %
    % function [fens,fes] = merge_nodes(fens, fes, tolerance)
    %
    % Merge by gluing together nodes from a single node set located within
    % tolerance of each other. The nodes are glued together by merging the
    % nodes that fall within a box of size "tolerance". The merged node
    % set, fens, and the finite element set with renumbered  connectivities
    % are returned.
    %
    % See also: fuse_nodes
    
    xyz1 = fens.xyz;
    dim =size(xyz1,2); 
    id1 = (1:count(fens))';
    c1=ones(size(xyz1,1),1);
    % Mark nodes from the array that are duplicated 
    for i=1:count(fens)
        if (id1(i)>0)% This note has not yet been marked for merging
            XYZ =xyz1(i,:);
            xyzd=abs(xyz1-c1*XYZ);% find the distances along  coordinate directions
            jx=find(sum(xyzd')'<tolerance);
            if (~isempty(jx))
                minn=min(jx);
                id1(jx) =-minn;
                id1(minn) =minn;
            end
        end
    end
    % Generate  merged arrays of the nodes
    xyzm = zeros(count(fens),dim);
    idm = zeros(count(fens),1);
    mid=1;
    for i=1:count(fens) % and then we pick only non-duplicated fens1
        if id1(i)>0
            id1(i)=mid;
            idm(mid)=mid;
            xyzm(mid,:)=xyz1(i,:);
            mid=mid+1;
        else 
            id1(i)=id1(-id1(i));
        end
    end
    nfens =mid-1;
    xyzm =xyzm(1:nfens,:);
    idm =idm(1:nfens,:);
    % Renumber the cells
    conns=fes.conn;
    for i=1:count(fes)
        conn=conns(i,:);
        conns(i,:)=id1(conn);
    end
    fes.conn=conns;
    clear fens
    fens=fenode_set(struct('xyz',xyzm(1:nfens,:)));
end

% function [fens,fes] = merge_nodes(fens, fes, tolerance)
%     % Merge together  nodes of a single node set.
%     %
%     % function [fens,fes] = merge_nodes(fens, fes, tolerance)
%     %
%     % Merge by gluing together nodes from a single node set located within
%     % tolerance of each other. The nodes are glued together by merging the
%     % nodes that fall within a box of size "tolerance". The merged node
%     % set, fens, and the finite element set with renumbered  connectivities
%     % are returned.
%     %
%     % See also: fuse_nodes
%
%     xyz1 = fens.xyz;
%     dim =size(xyz1,2);
%     % Mark nodes from the array that are duplicated
%     id1 = (1:count(fens))';
%     mgl=helper_mergegroups(fens.xyz, tolerance);
%     for c=1:length(mgl)
%         minnn=min(mgl{c}.n);
%         for k=1:length(mgl{c}.n)
%             if (mgl{c}.n(k)~=minnn)
%                 id1(mgl{c}.n(k)) =-minnn;
%             end
%         end
%     end
%     % Generate  merged arrays of the nodes
%     xyzm = zeros(count(fens),dim);
%     idm = zeros(count(fens),1);
%     mid=1;
%     for i=1:count(fens) % and then we pick only non-duplicated fens1
%         if id1(i)>0
%             id1(i)=mid;
%             idm(mid)=mid;
%             xyzm(mid,:)=xyz1(i,:);
%             mid=mid+1;
%         else
%             id1(i)=id1(-id1(i));
%         end
%     end
%     nfens =mid-1;
%     xyzm =xyzm(1:nfens,:);
%     idm =idm(1:nfens,:);
%     % Renumber the cells
%     conns=fes.conn;
%     for i=1:count(fes)
%         conn=conns(i,:);
%         conns(i,:)=id1(conn);
%     end
%     fes.conn=conns;
%     clear fens
%     fens=fenode_set(struct('xyz',xyzm(1:nfens,:)));
% end

