function [fens1,fes1] = mirror_mesh(fens, fes, Normal, Point, renumb)
% Mirror a 2-D mesh in a plane given by its normal and one point. 
%
% function [fens1,fes1] = mirror_mesh(fens, fes, Normal, Point, renumb)
%
% fens, gcells= mesh, 
% Normal, Point= 2-D arrays
%  
% Warning: The code to relies on the numbering of the cells: to reverse
% the orientation of the mirrored cells, the connectivity is listed in
% reverse order.   If the mirrored cells do not follow this rule (for instance 
% hexahedra for quadrilaterals), their areas will
% come out negative. In such a case the renumbering function 
% of the connectivity needs to be supplied..
%
% For instance: H8 elements require  the renumbering function to be supplied as
% [fens1,gcells1] = mirror_mesh(fens, gcells,...
%           [-1,0,0], [0,0,0], @(c)c([1, 4, 3, 2, 5, 8, 7, 6]));
% 
% 
    Normal = Normal/norm (Normal);
    fens1 =fens;
    xyz =fens1.xyz;
    for i=1:count(fens1)
        xyz(i,:) =xyz(i,:)-2*dot(Normal,xyz(i,:))*Normal;
    end
    fens1.xyz=xyz;
    % Reconnect the cells
    if (~exist('renumb','var'))
        renumb=@(conn)conn(end:-1:1);
    end
    fes1=fes;
    conns=fes1.conn;
    for i=1:size(conns,1)
        conns(i,:)=renumb(conns(i,:));
    end
    fes1.conn=conns;
end
