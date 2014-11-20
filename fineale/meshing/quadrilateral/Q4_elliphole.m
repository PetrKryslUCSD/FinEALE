function [fens,fes]=Q4_elliphole(xradius,yradius,L,H,nL,nH,nR,options)
% Mesh of one quarter of a rectangular plate with an elliptical hole
%
% function [fens,fes]=Q4_elliphole(xradius,yradius,L,H,nL,nH,nR,options)
%
% xradius,yradius = radius of the ellipse,
% L,H= and dimensions of the plate,
% nL,nH= numbers of edges along the side of the plate,
% nR= number of edges along the circumference,
% options= options accepted by fe_set_Q4
%
% Examples: 
%     [fens,fes]=Q4_elliphole(1.2,2.4,4.8,3.5,4,2,2,[]);
%     drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
%
%     [fens,fes]=Q4_elliphole(2.4,1.2,4.8,3.5,4,2,6,[]);
%     drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on

dA =pi/2/(nL +nH);
    tolerance =(xradius+yradius)/(nL*nH)/100;
    fens= []; fes= [];
    for i= 1:nH
        xy = [xradius*cos((i-1)*dA),yradius*sin((i-1)*dA);...
            L,(i-1)/nH*H;...
            L,(i)/nH*H;...
            xradius*cos((i)*dA),yradius*sin((i)*dA)];
        [fens1,fes1] = Q4_quadrilateral(xy,nR,1,options);
        if isempty(fens)
            fens=fens1; fes =fes1;
        else
            [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, tolerance);
            fes =cat(fes1,fes2);
        end
    end
    for i= 1:nL
        xy = [xradius*cos((nH+i-1)*dA),yradius*sin((nH+i-1)*dA);...
            (nL-i+1)/nL*L,H;...
            (nL-i)/nL*L,H;...
            xradius*cos((nH+i)*dA),yradius*sin((nH+i)*dA)];
        [fens1,fes1] = Q4_quadrilateral(xy,nR,1,options);
        [fens,fes1,fes2] = merge_meshes(fens1, fes1, fens, fes, tolerance);
        fes =cat(fes1,fes2);
    end
end
