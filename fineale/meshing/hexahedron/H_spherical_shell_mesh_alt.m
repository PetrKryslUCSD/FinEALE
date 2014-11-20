function [fens,fes]=H_spherical_shell_mesh_alt(radius,thickness,nrefine,nlayers, Convert)
% Create a hexahedral mesh of 1/8 of spherical shell.
%
% Create a volume mesh of 1/8 of the spherical shell of "radius" and "thickness". 
%
% function [fens,fes]=H_spherical_shell_mesh_alt(radius,thickness,nrefine,nlayers, Convert)
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
%     [fens,fes]=H_spherical_shell_mesh_alt(1.35,0.02,1,1);
%     drawmesh({fens,fes},'nodes','fes','facecolor','none'); hold on
%
%     [fens,fes]=H_spherical_shell_mesh_alt(1.35,0.02,1,1,@H8_to_H64);
%     drawmesh({fens,fes},'nodes','fes','facecolor','none'); hold on
%
% See also: 
%

[fens,fes]= H8_block(thickness, 1.0, 1.0, nlayers, 2*nrefine, 2*nrefine);
if exist('Convert','var') && ~isempty(Convert)
    [fens,fes] = Convert (fens,fes);
end
X=fens.xyz;
for j=1:count(fens)
    xp=[(radius+X(j,1))*cos(pi/2*X(j,2)),0,(radius+X(j,1))*sin(pi/2*X(j,2))];
    X(j,:) =[xp(1)*cos(pi/2*(1-X(j,3))),xp(1)*sin(pi/2*(1-X(j,3))),xp(3)];
end
fens.xyz=X;
tolerance =min([thickness/nlayers,radius/ (2*nrefine)]) /100;
[fens,fes] = merge_nodes(fens, fes, tolerance);
end
