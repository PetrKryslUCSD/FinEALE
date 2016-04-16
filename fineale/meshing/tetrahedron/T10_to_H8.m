function [fens,fes] = T10_to_H8(fens,fes)
% Convert a mesh of tetrahedra T10 (ten-node) to H8 hexahedra.
%
% function [fens,fes] = T10_to_H8(fens,fes)
%
% Each tetrahedron is split into four hexahedra, mid-face and interior
% nodes are inserted.
%
% Examples: 
% [fens,fes] = T4_sphere(3.1,1);
% [fens,fes] = T4_to_T10(fens,fes);
% [fens,fes] = T10_to_H8(fens,fes);
% fens= onto_sphere(fens,3.1,connected_nodes(mesh_boundary(fes,[])));
% figure; drawmesh({fens,fes},'fes','facecolor','y'); hold on
%
    nedges=6;
    emn = [5; 6; 7; 8; 9; 10];
    nfaces =4;
    fc = [1, 3, 2;1, 2, 4;2, 3, 4;1, 4, 3];
    Volumes = [1,5,11,7,8,12,15,14;
        5,2,6,11,12,9,13,15;
        12,9,13,15,8,4,10,14;
        7,11,6,3,14,15,13,10];
    nfe=[5,6,7; 5,9,8; 6,10,9; 7,10,8];
    conns = fes.conn; label=fes.label;
    nlabel=[];
    if (~isempty(label))
        nlabel=reshape(repmat(label',4,1),[],1);
    end
    n=count(fens);
    txyz=fens.xyz;
    % Preallocate the node arrays: the size is a guess
    EstimatedAddedn=ceil(size(conns,1)*11);
    txyz=[txyz;zeros(EstimatedAddedn,3)];
    cn=n;% Current node number
    clear ev i J other anchor exyz;
    % Make nodes at the midpoints of faces
    faces={};
    for i= 1:size(conns,1)
        for J = 1:nfaces
            fv=conns(i,fc(J,:));
            anchormin=min(fv);        anchormax=max(fv);
            Third =setdiff(fv, [anchormin,anchormax]);
            if length(faces)<anchormin
                faces{anchormin}=[];
            end
            if isempty(Have_face(faces{anchormin}(:,:),anchormax, Third))
                fxyz =(txyz(fv(1),:)+txyz(fv(2),:)+txyz(fv(3),:))/3;
                cn=cn+1; txyz(cn,:)=fxyz;
                faces{anchormin}=[faces{anchormin}(:,:); [anchormax,Third,cn]];
            end
        end
    end
    function ix2= Have_face(A, anchormax2, Third2)
        for j1=1:size(A,1)
            if (A(j1,1)==anchormax2)&&(A(j1,2)==Third2)
                ix2= j1; return;
            end
        end
        ix2= [];
    end
    function fen=find_face_node(fv3)
        anchormin3=min(fv3);        anchormax3=max(fv3);
        Third3 =setdiff(fv3, [anchormin3,anchormax3]);
        if isempty(faces{anchormin3})
            error(' No such face');
        end
        ix3=Have_face(faces{anchormin3}(:,:),anchormax3, Third3);
        if isempty(ix3)
            error(' No such face');
        end
        fen= faces{anchormin3}(ix3,3);
    end
    clear fv i J anchormin anchormax Third fxyz;
    % Add interior (volume) nodes
    interiors=zeros(size(conns,1),1);
    for i= 1:size(conns,1)
        vv=conns(i,:);
        vxyz =(txyz(vv(1),:)+txyz(vv(2),:)+txyz(vv(3),:)+txyz(vv(4),:))/4;
        cn=cn+1; txyz(cn,:)=vxyz;
        interiors(i)=cn;
    end
    function fen= find_volume_node(i4)
        fen= interiors(i4);
    end
    
    % construct new geometry cells
    nconns=zeros(4*count(fes),8);
    nc=1;
    for i= 1:size(conns,1)
        conn = [conns(i,:),zeros(1,5)];
        ix=5;
        for J = 1:nedges
            conn(ix)= conn(emn(J));
            ix=ix+1;
        end
        for J = 1:nfaces
            fv=conn(fc(J,:));
            conn(ix)=find_face_node(fv);
            %             Compute the position of the face node from the three edge nodes
            vxyz =mean(txyz(conn(nfe(J,:)),:));
            txyz(conn(ix),:)=vxyz;
            ix=ix+1;
        end
        conn(ix)=find_volume_node(i);
        % Compute the position of the interior node from the six edge nodes
        vxyz =mean(txyz(conn(end-6:end-1),:));
        txyz(conn(end),:)=vxyz;
        for J=1:4
            nconns(nc,:) =conn(Volumes(J,:));
            nc= nc+ 1;
        end
    end
    %Create new data structures
    fens=fenode_set(struct('xyz',txyz(1:cn,:)));
    fes =fe_set_H8(struct ('conn',nconns));
    if (~isempty(label))
        fes.label=nlabel;
    end
    return;
end

