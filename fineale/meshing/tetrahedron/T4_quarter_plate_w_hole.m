function [fens,fes] = T4_quarter_plate_w_hole(InternalRadius, Width, Thickness, nR, nW, nL, blockh)
% T4 mesh of a solid hollow cylinder.
% 
% [fens,fes] = T4_quarter_plate_w_hole(InternalRadius, Width, Length, nC, nT, nL, blockh)
%     
% InternalRadius= radius of the circular hole,
% Width = width of the plate
% Length=  length of the cylinder
% nR, nW, nL=nnumber of elements radially, through the
%      thickness, and  along the length
% blockh  = handle to a block-generating function
%
% Output:
% fens= finite element node set
% fes = finite element set
%
% Examples: 
%     [fens,fes] = T4_quarter_plate_w_hole(1.0, 2.5, 1.5, 7, 2, 3)
%     drawmesh({fens,fes},'fes','facecolor','red'); hold on
    
if (~exist('blockh','var'))
    blockh=@T4_blocka;
end
    [fens,fes] = blockh(2*pi/4, Width-InternalRadius, Thickness, 2*nW, nR, nL);
    fens = transform_apply(fens,@(x, data) (x+ [0, InternalRadius, 0]), []);
    climbPerRevolution= 0;
    fens = transform_2_helix(fens,climbPerRevolution);
    fens = rotate_mesh(fens,[0,0,pi/2],[0,0,0]);
    for  j=1:count(fens)
        d=fens.xyz(j,1:2);
        dn2=norm(d,2);
        dinf=norm(d,inf);
        r=(dn2-InternalRadius)/(Width-InternalRadius);
        fens.xyz(j,1:2)  =fens.xyz(j,1:2)+(-0.5*r+3/2*r^2)*(Width-dinf)*d/dinf;
    end    
end