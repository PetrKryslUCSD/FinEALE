function [fens] = translate_mesh(fens, TranslationVector)
% Translate a mesh along a vector. 
%
% function [fens] = translate_mesh(fens, TranslationVector)
%
% fens= Finite element nodes of the mesh, 
% TranslationVector= 3-D vector
% 
% 
    xyz =fens.xyz;
    shift=ones(size(xyz, 1),1)*reshape(TranslationVector,1,3);
    xyz = xyz+shift;
    fens.xyz=xyz;
end
