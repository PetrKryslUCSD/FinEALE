function [fens,fes] = L2_to_L3(fens,fes,options)
%Convert a mesh of L2 elements (two-node) to L3 line elements.
%
% function [fens,fes] = L2_to_L3(fens,fes,options)
%
% options =attributes recognized by the constructor of L3
%
% Examples: 
%     [fens,fes] = L2_block(5.0,4,struct('other_dimension',.25));
%     [fens,fes] = L2_to_L3(fens,fes,struct('other_dimension',.25))
%     drawmesh({fens,fes},'nodes','fes','facecolor','none', 'linewidth',3); hold on
%
% See also: L2_blockx, fe_set_L3


    if ~isstruct(options)
        options = struct('conn',[]);
    end
    % now generate new node number for each Element
    n=count(fens);
    ng =count(fes);
    id=(1:count(fens))';
    x =fens.xyz;
    id = [id;zeros(ng,1)];
    x = [x;zeros(ng,size(x,2))];
    conns =zeros(ng,3);
    conns(:,1:2) = fes.conn;
    conns(:,3)=n+(1:ng)';
    for i= 1:ng
        id(n+i)=n+i;
        x(n+i,:) =(x(conns(i,1),:)+x(conns(i,2),:))/2;
    end
    options.id =1:ng;
    options.conn =conns;
    options.other_dimension =fes.other_dimension;
    fes=fe_set_L3(options);
    fens=fenode_set(struct('xyz',x));
end
