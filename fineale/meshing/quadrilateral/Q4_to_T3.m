% Convert a mesh of quadrilateral Q4's to two T3 triangles  each.
%
% function [fens,fes] = Q4_to_T3(fens,fes,options)
%
% options =attributes recognized by the constructor fe_T3, and
% orientation = 'default' or 'alternate' chooses which diagonal is taken
%      for splitting
% Example: 
%     [fens,fes] = Q4_quadrilateral([-1,-1;2,-2;3,3;-1,1],2,3,[]);
%     [fens,fes] = Q4_to_T3(fens,fes,[]);
%     drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
%
%     [fens,fes] = Q4_quadrilateral([-1,-1;2,-2;3,3;-1,1],2,3,[]);
%     [fens,fes] = Q4_to_T3(fens,fes,struct('orientation','alternate'));
%     drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
% will generate triangles by using the alternate diagonal for splitting.
% 
% See also: Q4_to_T3_sd

function [fens,fes] = Q4_to_T3(fens,fes,options)
    if ~isstruct(options)
        options = struct('conn',[]);
    end
    connl1=(1:3);
    connl2=[1, 3, 4];
    if isfield(options,'orientation')
         if strcmp(options.orientation,'alternate')
             connl1=[1, 2, 4];
             connl2=[3, 4, 2];
         end
    end
    nedges=4;
    nconns=zeros(2*count(fes),3);
    conns = fes.conn;
    nc=1;
    for i= 1:count(fes)
        conn = conns(i,:);
        nconns(nc,:) =conn(connl1);
        nc= nc+ 1;
        nconns(nc,:) =conn(connl2);
        nc= nc+ 1;
    end
    options.conn=nconns;
    fes=fe_set_T3(options);
end
