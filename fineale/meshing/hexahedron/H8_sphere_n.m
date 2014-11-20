function [fens,fes]=H8_sphere_n(radius,nperradius)
% Create a solid mesh of 1/8 of sphere. 
%
% Create a solid mesh of 1/8 of the sphere of "radius", with nperradius 
% elements per radius. 
%
% function [fens,fes]=H8_sphere_n(radius,nperradius)
%
% Create a mesh of 1/8 of the sphere of "radius". The  mesh will consist of
% 4*(nperradius/2)^2 hexahedral elements.
% 
% Output:
% fens= finite element node set
% fes = finite element set
%
% Examples: 
%     [fens,fes]=H8_sphere_n(22.3,3);
%     drawmesh({fens,fes},'fes','facecolor','red'); hold on
%
% See also: H8_sphere


    if (mod(nperradius,2)~=0)
        nperradius=nperradius+1;
    end
    nL=nperradius/2; nW=nperradius/2; nH=nperradius/2;
    a=sqrt(2)/2;
    b=1/sqrt(3);
    c=0.6*a;
    d=0.6*b;
    xyz= [0,0,0;
        0.5,0,0;
        c,c,0;
        0,0.5,0;
        0,0,0.5;
        c,0,c;
        d,d,d;
        0,c,c;
        0,0,1;
        a,0,a;
        1,0,0;
        a,a,0;
        0,1,0;
        0,a,a;
        b,b,b]*radius;
    conn=[1,2,3,4,5,6,7,8;
        2,11,12,3,6,10,15,7;
        4,3,12,13,8,7,15,14;
        5,6,7,8,9,10,15,14];
    [fens,fes] = H8_hexahedron(xyz(conn(1,:),:),nL,nW,nH);
    [fens1,fes1] = H8_hexahedron(xyz(conn(2,:),:),nL,nW,nH);
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, norm(xyz)/1000);
    fes=cat(fes1,fes2);
    [fens1,fes1] = H8_hexahedron(xyz(conn(3,:),:),nL,nW,nH);
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, norm(xyz)/1000);
    fes=cat(fes1,fes2);
    [fens1,fes1] = H8_hexahedron(xyz(conn(4,:),:),nL,nW,nH);
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, norm(xyz)/1000);
    fes=cat(fes1,fes2);
    
    %     bg=mesh_boundary(fes);
    %     l=fe_select(fens,bg,struct ('facing', true, 'direction', [1,1,1]));
    %     cn = connected_nodes(subset(bg,l))   ;
    %     xyz=get(fens,'xyz');
    %     for j=1:length(cn)
    %         xyz(cn(j),:)=xyz(cn(j),:)*radius/norm(xyz(cn(j),:));
    %     end
    %     fens=set(fens,'xyz',xyz); 
    
    xyz=fens.xyz;
    layer=inf+zeros(size(xyz, 1),1);
    conn=fes.conn;
    bg=mesh_boundary(fes);
    l=fe_select(fens,bg,struct ('facing', true, 'direction', [1,1,1]));
    cn = connected_nodes(subset(bg,l))   ;
    layer(cn)=1;
    for j=1:nperradius-1
        for k=1:size(conn,1)
            ll=layer(conn(k,:));
            ml=min(ll);
            if (ml==j)
                ix=isinf(ll);
                ll(ix)=j+1;
                layer(conn(k,:))=ll;
            end
        end 
    end
    nxyz=xyz;
    for j=1:size(xyz,1)
        if (~isinf(layer(j)))
            nxyz(j,:)=nxyz(j,:)*(nperradius-layer(j)+1)/nperradius*radius/norm(nxyz(j,:));
        end
    end
    s= linspace(0, 1, nperradius);
    for j=1:size(xyz,1)
        if (~isinf(layer(j)))
            nxyz(j,:)=s(layer(j))*xyz(j,:)+(1-s(layer(j)))*nxyz(j,:);
        end
    end
    fens.xyz=nxyz;
end
