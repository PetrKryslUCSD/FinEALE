%Convert a mesh of L2 elements (two-node) to L4 line elements.
%
% function [fens,fes] = L2_to_L4(fens,fes,options)
%
% options =attributes recognized by the constructor of L3
%
% Examples: 
%     [fens,fes] = L2_block(5.0,4,struct('other_dimension',.25));
%     [fens,fes] = L2_to_L4(fens,fes,struct('other_dimension',.25))
%     drawmesh({fens,fes},'nodes','fes','facecolor','none', 'linewidth',3); hold on
%
% See also: L2_blockx, fe_set_L4
function [fens,fes] = L2_to_L4(fens,fes,options)
    if ~isstruct(options)
        options = struct('conn',[]);
    end
    % now generate new node number for each Element
    n=count(fens);
    ng =count(fes);
    x =fens.xyz;;
    x = [x;zeros(2*ng,size(x,2))];
    conns =zeros(ng,4);
    conns(:,[1,4]) = fes.conn;
    na=1;
    for i= 1:ng
        id(n+na)=n+na;
        x(n+na,:) =(2*x(conns(i,1),:)+x(conns(i,4),:))/3;
        conns(i,2)= n+na;
        na=na+1;
        x(n+na,:) =(x(conns(i,1),:)+2*x(conns(i,4),:))/3;
        conns(i,3)= n+na;
        na=na+1;
    end
    options.conn =conns;
    options.other_dimension = fes.other_dimension;
    fes=fe_set_L4(options);
    fens=fenode_set(struct('xyz',x));
end
