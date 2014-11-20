function [fens,fes] = Q4_refine_5(fens,fes)
% Refine a mesh of quadrilaterals into five interior quadrilaterals.
%
% function [fens,fes] = Q4_refine_5(fens,fes)
%
% Examples: 
% [fens,fes] = Q4_quadrilateral([-1,-1;2,-2;3,3;-1,1],2,3,[]);
% [fens,fes] = Q4_refine_5(fens,fes);
% drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
%
% See also: Q4_refine


    nsubcells=5;
    conns=fes.conn;
    xyz =fens.xyz;
    nfens=size(xyz,1);
    nifens=4*count(fes);% number of internal nodes
    xyz=[xyz;zeros(nifens,size(xyz,2))];
    % construct new geometry cells
    nconns =zeros(nsubcells*count(fes),4);
    nc=1;
    for i= 1:count(fes)
        xyzm=mean(xyz(conns(i,:),:));
        inn=nfens+(i-1)*4+(1:4);
        xyz(inn(1),:)=mean([xyzm;xyz(conns(i,1),:)]);
        xyz(inn(2),:)=mean([xyzm;xyz(conns(i,2),:)]);
        xyz(inn(3),:)=mean([xyzm;xyz(conns(i,3),:)]);
        xyz(inn(4),:)=mean([xyzm;xyz(conns(i,4),:)]);
        nconns(nc,:) =[conns(i,1) conns(i,2) inn(2) inn(1)];
        nc= nc+ 1;
        nconns(nc,:) =[conns(i,2) conns(i,3) inn(3) inn(2)];
        nc= nc+ 1;
        nconns(nc,:) =[conns(i,3) conns(i,4) inn(4) inn(3)];
        nc= nc+ 1;
        nconns(nc,:) =[conns(i,4) conns(i,1) inn(1) inn(4)];
        nc= nc+ 1;
        nconns(nc,:) =inn;
        nc= nc+ 1;
    end
    Thickness =fes.other_dimension;
    fes =fe_set_Q4(struct('conn',nconns,...
        'other_dimension',Thickness));
    fens =fenode_set(struct('xyz',xyz));
end
