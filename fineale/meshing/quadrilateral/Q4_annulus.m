function [fens,fes] = Q4_annulus(rin,rex,nr,nc,thickness)
% Mesh of an annulus segment.
%
% function [fens,fes] = Q4_annulus(rin,rex,nr,nc,thickness)
%
% Mesh of an annulus segment in the first quadrant, centered at the origin,
% with internal radius rin, and  external radius rex. Divided into
% elements: nr, nc in the radial and circumferential direction
% respectively.
%
% Examples: 
%     [fens,fes] = Q4_annulus(1.5,2.5,2,17,1.0);
%     drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
%
% See also: Q4_block
%
trin=min(rin,rex);
trex=max(rin,rex);
[fens,fes]=Q4_block(trex-trin,pi/2,nr,nc, thickness);
xy=fens.xyz;
for i=1:count (fens)
    r=trin+xy(i,1); a=xy(i,2);
    xy(i,:)=[r*cos(a) r*sin(a)];
end
fens.xyz=xy;
end
