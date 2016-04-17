% Test generation of a hollow cylinder And conversion of quadratic
% tetrahedra  into linear Hexahedra.
function test_hcyl_T10_to_H8
Thickness= 0.2;
Length= 2.5;
InternalRadius=3* Thickness;
nC=4; nT=1; nL=3;
[fens,fes] = T10_block(pi,Thickness, Length, nC, nT, nL);
fens = transform_apply(fens,@(x, data) (x+ [0, (cos(2*x(3))+2)/2*InternalRadius, 0]), []);
climbPerRevolution= 0;
fens = transform_2_helix(fens,climbPerRevolution);
[fens,fes] = T10_to_H8(fens,fes);
bg=mesh_boundary(fes);
gv =drawmesh({fens,bg},'fes','facecolor','y');
headlight(gv);

    v= hex_volumes(fens.xyz, fes);
    
    any(v< 0.0)
    
function v= hex_volumes(xs, gcells);
    integration_rule = gauss_rule(struct('dim',3,'order' ,2));
    pc = integration_rule.param_coords;
    w  = integration_rule.weights;
    
    npts_per_gcell = integration_rule.npts;
    conns = fes.conn; % connectivity
    v=zeros(size(conns,1),1);
    for i=1:size(conns,1)
        conn =conns(i,:);
        x=xs(conn,:);
        for j=1:npts_per_gcell
            N = bfun(gcells,pc(j,:));
            Nder = bfundpar(gcells,pc(j,:));
            J = Jacobian_matrix(gcells,Nder,x);
            Jac = Jacobian_volume(gcells,conn, N, J, x); 
            v(i)=v(i)+Jac;
        end
    end
end

end