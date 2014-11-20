% Surface mesh of 1/8 of a sphere with a given number of elements per circumference. 
%
% function [fens,fes]=Q4_sphere_n(radius,nperradius,thickness)
%
% Create a mesh of 1/8 of the sphere of "radius". The  mesh will consist of
% 3*(nperradius/2)^2 quadrilateral elements, which corresponds to 4*nperradius
% element edges per circumference. The parameter nperradius should be an even 
% number; if that isn't so is adjusted to by adding one. 
% 
% The quadrilaterals have thickness "thickness".
%
% Examples: 
% [fens,fes]=Q4_sphere_n(77.1,5,1.0);
% drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
%
% See also: Q4_sphere

function [fens,fes]=Q4_sphere_n(radius,nperradius,thickness)
    if (mod(nperradius,2)~=0)
        nperradius=nperradius+1;
    end
    nL=nperradius/2; nW=nperradius/2; 
    a=sqrt(2)/2;
    b=1/sqrt(3);
    c=0.6*a;
    d=0.6*b;
    xyz=[[1, 0, 0];[0, 1, 0];[0, 0, 1];[a, a, 0];[0, a, a];[a, 0, a];[b, b, b]];
    conn=[[1, 4, 7, 6];[4, 2, 5, 7]; [3, 6, 7, 5]];
    [fens,fes] = Q4_quadrilateral(xyz(conn(1,:),:),nL,nW,struct ('other_dimension', thickness));
    [fens1,fes1] = Q4_quadrilateral(xyz(conn(2,:),:),nL,nW,struct ('other_dimension', thickness));
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, norm(xyz)/1000);
    fes=cat(fes1,fes2);
    [fens1,fes1] = Q4_quadrilateral(xyz(conn(3,:),:),nL,nW,struct ('other_dimension', thickness));
    [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, norm(xyz)/1000);
    fes=cat(fes1,fes2);
    fens= onto_sphere(fens, radius);
end
