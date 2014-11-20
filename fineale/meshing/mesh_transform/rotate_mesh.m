function [fens] = rotate_mesh(fens, RotationVector, Point)
% Rotate a mesh around a vector anchored at the given point. 
%
% function [fens] = rotate_mesh(fens, RotationVector, Point)
%
% fens= Finite element nodes of the mesh, 
% RotationVector, Point= 3-D vectors
% 
% 
    R = rotmat(RotationVector);
    xyz =fens.xyz;
    pivot=ones(size(xyz, 1),1)*reshape(Point,1,3);
    xyz = pivot+(R*(xyz-pivot)')';
    fens.xyz=xyz;
end
