function [fens,fes] = H8_hollow_cylinder(InternalRadius, Thickness, Length, nC, nT, nL)
% H8 mesh of a solid hollow cylinder.
% 
% function [fens,fes] = H8_hollow_cylinder(InternalRadius, Thickness, Length, nC, nT, nL)
%     
% Output:
% fens= finite element node set
% fes = finite element set
%
% Examples: 
%     [fens,fes] = H8_hollow_cylinder(1.0, 00.5, 1.5, 17, 2, 3)
%     drawmesh({fens,fes},'fes','facecolor','red'); hold on
    
    [fens,fes] = H8_block(2*pi,Thickness, Length, nC, nT, nL);
    fens = transform_apply(fens,@(x, data) (x+ [0, InternalRadius, 0]), []);
    climbPerRevolution= 0;
    fens = transform_2_helix(fens,climbPerRevolution);
    [fens,fes] = merge_nodes(fens, fes, min( [Thickness/nT,Length/nL,2*pi*InternalRadius/nC] )/100);
end