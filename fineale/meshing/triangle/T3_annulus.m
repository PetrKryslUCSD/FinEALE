function [fens,fes] = T3_annulus(rin,rex,nr,nc,thickness)
% Mesh of an annulus segment, triangles T3.
%
% function [fens,fes] = T3_annulus(rin,rex,nr,nc,thickness)
%
% Mesh of an annulus segment in the first quadrant, centered at the origin,
% with internal radius rin, and  external radius rex. 
% Divided into elements: nr, nc in the radial and circumferential direction respectively.
% options = Attributes recognized by the constructor of fe_set_T3.
%
% Examples: 
%     [fens,fes] = T3_annulus(11.1,42.2,3,12,[]);
%     figure; drawmesh({fens,fes},'fes','facecolor','y'); hold on
%
% See also:  T3_blocku, T3_cblock,  T3_cblock2d, T3_crossblock, T3_ablock,
%            T3_ablockc, T3_ablockn, T3_block, T3_blockc, T3_blockn

trin=min(rin,rex);
trex=max(rin,rex);
[fens,fes]=T3_block(trex-trin,pi/2,nr,nc, thickness);
xy=fens.xyz;
for i=1:count(fens)
    r=trin+xy(i,1); a=xy(i,2);
    xy(i,:)=[r*cos(a) r*sin(a)];
end
fens.xyz=xy;