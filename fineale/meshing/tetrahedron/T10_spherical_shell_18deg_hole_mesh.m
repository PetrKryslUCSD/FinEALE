function [fens,fes]=T10_spherical_shell_18deg_hole_mesh(radius,thickness,ncirc,nmerid,nlayers, blockh)
% Create a tetrahedral mesh of 1/8 of spherical shell.
%
% Create a volume mesh of 1/8 of the spherical shell of "radius" and "thickness". 
%
% function [fens,fes]=H_spherical_shell_mesh(radius,thickness,nrefine,nlayers, Convert)
%
% Create a hexahedral mesh of 1/8 of the spherical shell of "radius"  and
% "thickness". The radius is the _internal_ radius.  The external radius is
% radius + thickness. The  mesh will consist of three*nlayers hexahedral
% elements if "nrefine==0", or more if "nrefine>0". "nrefine" is the number
% of bisections applied  to refine the mesh. The shell has the thickness
% "thickness".
%
% The mesh is first created to consist of H8 Hexahedra.  It may be then
% converted to other hexahedral types if Convert is a handle to a
% conversion function (for instance H8_to_H20).
%
% Examples:
%     [fens,fes]=T10_spherical_shell_18deg_hole_mesh(1.35,0.02,1,1);
%     drawmesh({fens,fes},'nodes','fes','facecolor','none'); hold on
%
%     [fens,fes]=T10_spherical_shell_18deg_hole_mesh(1.35,0.02,1,1,@H8_to_H64);
%     drawmesh({fens,fes},'nodes','fes','facecolor','none'); hold on
%
% See also: 
%
if (~exist('blockh','var'))
    blockh=@T10_blocka;
end
[fens,fes] =blockh(90,90-18,thickness,ncirc,nmerid,nlayers);
xy=fens.xyz;
for i=1:count (fens)
    p=xy(i,1)/180*pi; f=xy(i,2)/180*pi; r=radius+xy(i,3);
    xy(i,:)=[r*cos(p)*cos(f) r*sin(p)*cos(f) r*sin(f)];
end
fens.xyz=xy;
% drawmesh({fens,fes},'fes')
% disp()