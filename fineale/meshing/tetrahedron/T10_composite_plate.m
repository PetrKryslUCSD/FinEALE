function [fens,fes] = T10_composite_plate(L,W,ts,nL,nW,nts)
% T10 block mesh for a layered block (composite plate).
%
% function [fens,fes] = T10_composite_plate(L,W,ts,nL,nW,nts)
%
% L,W= length and width,
% ts= Array of layer thicknesses,
% nL,nW= Number of elements per length and width,
% nts= array of numbers of elements per layer
%
% The fes of each layer are labeled with the layer number.
%
% Output:
% fens= finite element node set
% fes = finite element set
%
%
% Examples: 
%     a=200; b=600; h=50;
%     angles =[0,90,0];
%     nLayers =length(angles);
%     na=4; nb=4;
%     nts= 1*ones(nLayers,1);% number of elements per layer
%     ts= h/nLayers*ones(nLayers,1);% layer thicknesses
%     [fens,fes] = H8_composite_plate(a,b,ts,na,nb,nts);;
%     gv=drawmesh( {fens,subset(fes,fe_select(fens,fes,struct('label', 1)))},'fes', 'facecolor','r');
%     gv=drawmesh( {fens,subset(fes,fe_select(fens,fes,struct('label', 2)))},'gv',gv,'fes', 'facecolor','g');
%     gv=drawmesh( {fens,subset(fes,fe_select(fens,fes,struct('label', 3)))},'gv',gv,'fes', 'facecolor','b');
%
%
% See also: T10_block
%
gtolerance =min(abs(ts))/max(nts)/10;
hblock=@T4_blocka;
totthel =sum(nts);
[fens,fes] = hblock(L,W,totthel*min(ts),nL,nW,totthel);
fes.label=0;
ofens=fens;
cnlayer =1; layerz=0; nthel=0;
for layer =1:length(ts)
    el=fe_select(ofens,fes,struct('box', [-inf,inf,-inf,inf,nthel*min(ts),(nthel+nts(layer))*min(ts)], 'inflate', gtolerance));
    label=fes.label; label(el)=layer;
    fes.label = label;;
    th=ts(layer);
    for slayer =1:nts(layer)
        nl=fenode_select(ofens,struct('box', [-inf,inf,-inf,inf,cnlayer*min(ts),cnlayer*min(ts)], 'inflate', gtolerance));
        fens.xyz(nl,3)=layerz+(slayer)*th/nts(layer);
        cnlayer=cnlayer+1;
    end    
    layerz=layerz+th;
    nthel=nthel+nts(layer);
end
 [fens,fes] = T4_to_T10(fens,fes);
end
