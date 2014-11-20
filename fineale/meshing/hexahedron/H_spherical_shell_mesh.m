function [fens,fes]=H_spherical_shell_mesh(radius,thickness,nrefine,nlayers, Convert)
% Create a hexahedral mesh of 1/8 of spherical shell.
%
% Create a volume mesh of 1/8 of the spherical shell of "radius" and "thickness". 
%
% function [fens,fes]=H_spherical_shell_mesh(radius,thickness,nrefine,nlayers, Convert)
%
% Create a hexahedral mesh of 1/8 of the spherical shell of "radius"  and
% "thickness". The radius is the _internal_ radius.  The external radius is
% radius + thickness. The  mesh will consist of three*nlayers hexahedral
% elements if "nrefine==0", or more if "nrefine>0". "nrefine" is the number
% of bisections applied  to refine the mesh. The shell has the thickness
% "thickness".
%
% The mesh is first created to consist of H8 Hexahedra.  It may be then
% converted to other hexahedral types if Convert is a handle to a
% conversion function (for instance H8_to_H20).
%
% Examples:
%     [fens,fes]=H_spherical_shell_mesh(1.35,0.02,1,1);
%     drawmesh({fens,fes},'nodes','fes','facecolor','none'); hold on
%
%     [fens,fes]=H_spherical_shell_mesh(1.35,0.02,1,1,@H8_to_H64);
%     drawmesh({fens,fes},'nodes','fes','facecolor','none'); hold on
%
% See also: 
%

p=1/3;
fens =fenode_set(struct('xyz',[0, 0; 1,0; 0,1; 1/2,0; 1/2,1/2; 0,1/2; 1/2*(1-p),p]));
    fes=fe_set_Q4(struct ('conn',[[1, 4, 7, 6];[4, 2, 5, 7]; [3, 6, 7, 5]],'other_dimension', thickness));
    for i=1:nrefine
        [fens,fes]=Q4_refine(fens,fes);
    end
    function xyz= z(xyz, layer)
        xyz= [xyz,layer/nlayers*thickness];
    end
	X=[+1,0,0; 0,+1,0; 0,0,+1; sqrt(2)/2,sqrt(2)/2,0; 0,sqrt(2)/2,sqrt(2)/2; sqrt(2)/2,0,sqrt(2)/2];
    [fens,fes] = H8_extrude_Q4(fens,fes,nlayers,@z);
	if exist('Convert','var') && ~isempty(Convert)
	    [fens,fes] = Convert (fens,fes);
	end
    xs= fens.xyz;
	Newxs=xs;
    for i=1:size(xs,1)
    	p= X'*T6bfun(xs(i,1:2));
		R= radius+xs(i,3);
    	Newxs(i,:)= R*p/ norm(p);
    end
    fens.xyz= Newxs;
end

function N = T6bfun(param_coords)
    r=param_coords(1);
    s=param_coords(2);
    t = 1 - r - s;
    N = [...
        t * (t + t - 1);
        r * (r + r - 1);
        s * (s + s - 1);
        4 * r * t;
        4 * r * s;
        4 * s * t];
end