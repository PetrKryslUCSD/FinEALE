function [fens,fes] = Q4_quadrilateral(xyz,nL,nW,options)
% Mesh of a general quadrilateral given by the location of the vertices.
%
% function [fens,fes] = Q4_quadrilateral(xyz,nL,nW,options)
%
% xyz = One vertex location per row; Either two rows (for a rectangular
% block given by the two corners), or four rows (General quadrilateral).
% Divided into elements: nL, nW in the first and second direction.
% options = Attributes recognized by the constructor of fe_set_Q4.
%
% Examples: 
% [fens,fes] = Q4_quadrilateral([-1,-1;2,-2;3,3;-1,1],2,3,[]);
% drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
%
    npts=size(xyz,1);
    if npts==2
        lo=min(xyz);
        hi=max(xyz);
        xyz=[[lo(1),lo(2)];...
            [hi(1),lo(2)];...
            [hi(1),hi(2)];...
            [lo(1),hi(2)]];
    elseif npts~=4
        error('Need 2 or 8 points');
    end


    [fens,fes] = Q4_block(2,2,nL,nW,options);

    xyz1=fens.xyz;
    if (size(xyz1,2)<size(xyz,2))
        nxyz1=zeros(size(xyz1,1),size(xyz,2));
        nxyz1(:,1:size(xyz1,2))=xyz1;
        xyz1=nxyz1;
    end
    for i=1:count(fens)
        N = bfun(xyz1(i,:)-1);% shift coordinates by -1
        xyz1(i,:) =N'*xyz;
    end
    fens.xyz=xyz1;
end

function val = bfun(param_coords)
    one_minus_xi = (1 - param_coords(1));
    one_plus_xi  = (1 + param_coords(1));
    one_minus_eta = (1 - param_coords(2));
    one_plus_eta  = (1 + param_coords(2));

    val = [0.25 * one_minus_xi * one_minus_eta;
        0.25 * one_plus_xi  * one_minus_eta;
        0.25 * one_plus_xi  * one_plus_eta;
        0.25 * one_minus_xi * one_plus_eta];
    return;
end
