% Draw a basis function on a 2D finite element mesh.
%
% function drawbasis2d(fens, fes, nodenum, scale)
%
%   where
%     fens     - array of nodes
%     fes   - array of geometry cells
%     nodenum  - number of the node whose basis function is to be
%                visualized
%     scale - how it should be scaled in height
%
function drawbasis2d(fens, fes, nodenum, scale)
% Mesh
if (get(fes(1),'dim') ~= 2)
    error(['Makes sense only for 2D meshes: got dim=' num2str(get(fes(1),'dim')) '!']);
    return;
end
if (nodenum < 1 | nodenum > length(fens))
    error(['Invalid node number: ' num2str(nodenum) ' (valid numbers are 1 < ... <= '...
        num2str(length(fens)) ')']);
    return;
end

% Geometry
geom = field(struct('name',['geom'], 'dim', 2, 'fens',fens));
xy=get(geom,'values');
geom = field(struct('name',['geom'], 'dim', 3, 'data', [xy 0*xy(:,1)] ));

u = 0 * clone(geom,'u');
u = scatter(u,[nodenum],[0 0 1]);
U = get (u,'values');

gv=graphic_viewer;
gv=reset (gv,struct('limits',[min(xy(:,1)) max(xy(:,1)) min(xy(:,2)) max(xy(:,2))  0 scale*1]));
dcm=data_colormap(struct ('range',[0, scale*1], 'colormap',jet));
colorfield=field(struct ('name', ['colorfield'], 'data',map_data(dcm, 0.8*scale*U(:, 3))));
for i=1:count(fes)
    draw(fes(i),gv, struct ('x', geom, 'u', 0*u, 'facecolor','none'));
    draw(fes(i),gv, struct ('x', geom, 'u', +scale*u,'colorfield',colorfield, 'shrink', 1.0));
end
return;
