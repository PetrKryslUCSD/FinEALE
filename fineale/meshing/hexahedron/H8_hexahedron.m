function [fens,fes] = H8_hexahedron(xyz,nL,nW,nH,block_mesh_handle)
% Mesh of a general hexahedron given by the location of the vertices.
%
% function [fens,fes] = H8_hexahedron(xyz,nL,nW,nH,block_mesh_handle)
% 
% xyz = One vertex location per row; Either two rows (for a rectangular
%      block given by the its corners), or eight rows (general hexahedron).
% nL, nW, nH = Divided into elements: nL, nW, nH in the first, second, and
%      third direction. 
% Optional argument:
% block_mesh_handle = function handle of the block-generating mesh function
%      (having the signature of the function H8_block()).
%
% Output:
% fens= finite element node set
% fes = finite element set
%
%
% Examples: 
%
%     xyz = [3, 1, 6; -5, 2, 1];
%     [fens,fes] = H8_hexahedron(xyz,12,3,4);
%     drawmesh({fens,fes},'fes','facecolor','red'); hold on
%
%     A=[0,0,0]; B=[0,0,2]; C=[0,3,2]; D=[0,3,0];
%     E=[5,0,0]; F=[5,0,2]; G=[5,3,2]; H=[5,3,0];
%     P=[3.75,0,0];
%     [fens,fes] = H8_hexahedron([A;P;(D+H)/2;D;B;(B+F)/2;(C+G)/2;C],2,3,4,[]);
%     drawmesh({fens,fes},'fes','facecolor','red'); hold on
%
%     A=[0,0,0]; B=[0,0,2]; C=[0,3,2]; D=[0,3,0];
%     E=[5,0,0]; F=[5,0,2]; G=[5,3,2]; H=[5,3,0];
%     P=[3.75,0,0];
%     [fens,fes] = H8_hexahedron([A;P;(D+H)/2;D;B;(B+F)/2;(C+G)/2;C],1,2,3,@H20_block);
%     drawmesh({fens,fes},'nodes','fes','facecolor','none'); hold on

    npts=size(xyz,1);
    if npts==2
        lo=min(xyz);
        hi=max(xyz);
        xyz=[[lo(1),lo(2),lo(3)];...
            [hi(1),lo(2),lo(3)];...
            [hi(1),hi(2),lo(3)];...
            [lo(1),hi(2),lo(3)];...
            [lo(1),lo(2),hi(3)];...
            [hi(1),lo(2),hi(3)];...
            [hi(1),hi(2),hi(3)];...
            [lo(1),hi(2),hi(3)]];
    elseif npts~=8
        error('Need 2 or 8 points');
    end

    if ~exist('block_mesh_handle') || isempty(block_mesh_handle)
        block_mesh_handle =@H8_block;
    end

    [fens,fes] = block_mesh_handle(2,2,2,nL,nW,nH);

    pxyz=fens.xyz;
    for i=1:count(fens)
        N = bfun(pxyz(i,:)-1);% shift coordinates by -1
        pxyz(i,:) =N'*xyz;
    end
    fens.xyz=pxyz;
end

function val = bfun(param_coords)
    one_minus_xi    = (1 - param_coords(1));
    one_minus_eta   = (1 - param_coords(2));
    one_minus_theta = (1 - param_coords(3));
    one_plus_xi     = (1 + param_coords(1));
    one_plus_eta    = (1 + param_coords(2));
    one_plus_theta  = (1 + param_coords(3));
    val = [one_minus_xi*one_minus_eta*one_minus_theta;...
        one_plus_xi*one_minus_eta*one_minus_theta;...
        one_plus_xi*one_plus_eta*one_minus_theta;...
        one_minus_xi*one_plus_eta*one_minus_theta;...
        one_minus_xi*one_minus_eta*one_plus_theta;...
        one_plus_xi*one_minus_eta*one_plus_theta;...
        one_plus_xi*one_plus_eta*one_plus_theta;...
        one_minus_xi*one_plus_eta*one_plus_theta] / 8;
    return;
end
