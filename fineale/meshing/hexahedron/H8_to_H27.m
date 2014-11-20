function [fens,fes] = H8_to_H27(fens,fes)
% Convert a mesh of hexahedra H8 to hexahedra H27.
%
% function [fens,fes] = H8_to_H27(fens,fes)
%
% Arguments and
% Output:
% fens= finite element node set
% fes = finite element set
%
    nedges=12;
    nfaces=6;
    ec = [1, 2; 2, 3; 3, 4; 4, 1; 5, 6; 6, 7; 7, 8; 8, 5; 1, 5; 2, 6; 3, 7; 4, 8;];
    fc = [1     4     3     2;
        1     2     6     5;
        2     3     7     6;
        3     4     8     7;
        4     1     5     8;
        6     7     8     5];
    conns = fes.conn;
    % make a search structure for edges
    edges={};
    for i= 1:count(fes)
        conn = conns(i,:);
        for J = 1:nedges
            ev=conn(ec(J,:));
            anchor=min(ev);
            if length(edges)<anchor
                edges{anchor}=[];
            end
            edges{anchor}=unique([edges{anchor} max(ev)]);
        end
    end
    % make a search structure for faces
    faces={};
    for i= 1:count(fes)
        conn = conns(i,:);
        for J = 1:nfaces
            fv=conn(fc(J,:));
            anchor=min(fv);
            if length(faces)<anchor
                faces{anchor}={};
            end
            sfv=sort(fv);
            sf=  faces{anchor};
            newf=[anchor sfv(2:end)];
            fnd=0;
            for kk=1:length(sf)
                if sf{kk}==newf
                    fnd=1; break;
                end
            end
            if ~fnd
                sf{end+1}=newf;
            end
            faces{anchor}=sf;
        end
    end
    %     number of vertex nodes
    nvn=count(fens);
    xyz1 =fens.xyz;
    id1= (1:size(xyz1,1))';
    %   now generate new node number for each edge
    n=count(fens);
    enodes=edges;
    for i= 1:length(edges)
        e=edges{i};
        for J = 1:length(e)
            n=n+1;
            e(J)=n;
        end
        enodes{i}=e;
    end%
    % Allocate for vertex nodes plus edge nodes
    xyz =zeros(n,3);
    xyz(1:size(xyz1,1),:) = xyz1;
    id =zeros(n,1);
    id(1:size(xyz1,1),:) = id1;
    % number of edge nodes
    nen=n-nvn;
    % calculate the locations of the new nodes
    % and construct the new nodes
    for i= 1:length(edges)
        e=edges{i};
        n=enodes{i};
        for J = 1:length(e)
            id(n(J)) =n(J);
            xyz(n(J),:)=mean(xyz([i,e(J)],:));
        end
    end
    % now generate new node number for each face
    n=size(xyz,1);
    fnodes=faces;
    for i= 1:length(faces)
        f=faces{i};
        fn=zeros(length(f),1);
        for J = 1:length(f)
            n=n+1;
            fn(J)=n;
        end
        fnodes{i}=fn;
    end
    % Preallocate to add face nodes
    xyz1 =xyz;
    id1 =id;
    xyz =zeros(n,3);
    xyz(1:size(xyz1,1),:) = xyz1;
    id =zeros(n,1);
    id(1:size(xyz1,1),:) = id1;
    % calculate the locations of the new nodes
    % and construct the new nodes
    for i= 1:length(faces)
        f=faces{i};
        fn=fnodes{i};
        for J = 1:length(f)
            nf=f{J};
            xyz(fn(J),:) = mean(xyz(nf(:),:));
            id(fn(J))=fn(J);
        end
    end
    % now generate new node number for each gcell
    vnodes=zeros(1,count(fes));
    n=size(xyz,1);
    for i= 1:count(fes)
        n=n+1;
        vnodes(i)=n;
    end
    % Preallocate to add face nodes
    xyz1 =xyz;
    id1 =id;
    xyz =zeros(n,3);
    xyz(1:size(xyz1,1),:) = xyz1;
    id =zeros(n,1);
    id(1:size(xyz1,1),:) = id1;
    % calculate the locations of the new nodes
    % and construct the new nodes
    for i= 1:count(fes)
        conn = conns(i,:);
        xyz(vnodes(i),:)=mean(xyz(conn,:));
        id(vnodes(i))=vnodes(i);
    end
    nfens=size(xyz,1);% number of nodes in the original mesh plus number of the edge nodes
    % construct new geometry cells
    nconns =zeros(size(conn,1),27);
    nc=1;
    for i= 1:count(fes)
        conn = conns(i,:);
        econn=zeros(1,nedges);
        for J = 1:nedges
            ev=conn(ec(J,:));
            anchor=min(ev);
            e=edges{anchor};
            n=enodes{anchor};
            econn(J)=n(find(e==max(ev)));
        end
        fconn=zeros(1,nfaces);
        for J = 1:nfaces
            fv=conn(fc(J,:));
            anchor=min(fv);
            f=faces{anchor};
            fn=fnodes{anchor};
            af=sort(fv);
            fnd=0;
            for kk=1:length(f)
                if f{kk}==af
                    fnd=1; break;
                end
            end
            if ~fnd
                error('Not found');
            end
            fconn(J)=fn(kk);
        end
        vconn=vnodes(i);
        nconns(nc,:) =[conn econn fconn vconn];
        nc= nc+ 1;
    end
    fens.xyz=xyz;
    fes=fe_set_H27(struct('conn',nconns)) ;
    return;
end
