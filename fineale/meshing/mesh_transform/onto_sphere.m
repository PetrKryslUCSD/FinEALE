function fens= onto_sphere(fens,radius,list)
% Project nodes onto a sphere of given radius.
%
% function fens= onto_sphere(fens,radius,list)
%
% fens = finite element node set,
% radius = radius of the sphere,
% list = optional argument, if not empty then only 
%	     the nodes in the list are to be moved; otherwise all nodes are moved.  
    xyz=fens.xyz;
    if (~exist( 'list', 'var' ))||isempty(list)
        list=1:size(xyz,1);
    end
    for i=1:  length(list)
        j= list(i);
        xyz(j,:) =xyz(j,:)*radius/norm(xyz(j,:));
    end
    fens.xyz=xyz;
end
