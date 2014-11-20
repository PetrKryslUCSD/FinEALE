function [fens,fes]=H8_sphere(radius,nrefine)
% Create a solid mesh of 1/8 of sphere. 
%
% Create a solid mesh of 1/8 of the sphere of "radius", refined nrefine-times. 
%
% function [fens,fes]=H8_sphere(radius,nrefine)
%
% Create a mesh of 1/8 of the sphere of "radius". The  mesh will consist of
% four hexahedral elements if "nrefine==0", or more if "nrefine>0".
% "nrefine" is the number of bisections applied  to refine the surface mesh.
% 
% Output:
% fens= finite element node set
% fes = finite element set
%
% Examples: 
%     [fens,fes]=H8_sphere(22.3,0);
%     drawmesh({fens,fes},'fes','facecolor','red'); hold on
%
%     [fens,fes]=H8_sphere(22.3,2);
%     drawmesh({fens,fes},'fes','facecolor','red'); hold on
%
% See also: H8_sphere_n

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
    fens =fenode_set(struct('xyz',xyz));
    fes=fe_set_H8(struct ('conn',conn));
    for i=1:nrefine
        [fens,fes]=H8_refine(fens,fes);
        bg=mesh_boundary(fes);
        l=fe_select(fens,bg,struct ('facing', true, 'direction', [1,1,1]));
        cn = connected_nodes(subset(bg,l))   ;
        xyz=fens.xyz;
        for j=1:length(cn)
            xyz(cn(j),:)=xyz(cn(j),:)*radius/norm(xyz(cn(j),:));
        end
        fens.xyz=xyz;
    end
end
