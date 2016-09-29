function [fens,fes] = T4_quarter_hollow_cylinder(InternalRadius, Thickness, Length, nC, nT, nL, blockh)
% T4 mesh of a solid hollow cylinder.
% 
% function [fens,fes] = T4_hollow_cylinder(InternalRadius, Thickness,
%          Length, nC, nT, nL, blockh)
%     
% InternalRadius, Thickness = internal radius and thickness of the wall
% Length=  length of the cylinder
% nC, nT, nL=nnumber of elements circumferentially, through the
%      thickness, and  along the length
% blockh  = handle to a block-generating function
%
% Output:
% fens= finite element node set
% fes = finite element set
%
% Examples: 
%     [fens,fes] = T4_quarter_hollow_cylinder(1.0, 00.5, 1.5, 17, 2, 3)
%     drawmesh({fens,fes},'fes','facecolor','red'); hold on
    
if (~exist('blockh','var'))
    blockh=@T4_blocka;
end

    [fens,fes] = blockh(2*pi/4,Thickness, Length, nC, nT, nL);
    fens = transform_apply(fens,@(x, data) (x+ [0, InternalRadius, 0]), []);
    climbPerRevolution= 0;
    fens = transform_2_helix(fens,climbPerRevolution);
end