function curved_beam_h20
%%
% Large-deflection problem. Solved many times:
%
% Bathe, Boulourchi 1979
% Simo, Vu-Quoc 1986
% Cardona, Geradin 1988
% Krysl 1993
% Krysl, FAESOR script curved_cantilever

% Present calculation refers to the data of
% Reference: EULERIAN FORMULATION FOR LARGE-DISPLACEMENT ANALYSIS OF SPACE
% FRAMES, by B.  A.  lzzuddin I and A.  S.  Elnashai,
% Journal of Engineering Mechanics, Vol.  119, No.  3,  March, 1993.

%%
%
% Parameters:
E=1e7;% lb.in^-2
nu=0.3;
h= 1; b =1;% in
Fmag=600;%  lb
radius=100;% in
ang=45;
nel=16;% number of elements
utol = h/1e9;
htol = h/100;
graphics = true;
nincr  =10;

% prop = property_deformation_neohookean (struct('E',E,'nu',nu));
% mater = material_deformation_neohookean_triax(struct('property',prop));
prop = property_deformation_linear_iso (struct('E',E,'nu',nu));
mater = material_deformation_stvk_triax(struct('property',prop));

surface_integration_rule=gauss_rule(struct('dim',2,'order',2));


%          Mesh
% Create the mesh and initialize the geometry
[fens,fes]= H20_block(h,ang/180*pi,b,2,nel,2);
bdry_fes = mesh_boundary(fes, []);
icl = fe_select(fens, bdry_fes, struct('box', [-100*radius, 100*radius,ang/180*pi,ang/180*pi,-100*radius, 100*radius],'inflate',htol));
xy=fens.xyz;
for i=1: size(xy,1)
    r=radius-h/2+xy(i,1); a=xy(i,2);
    xy(i,:)=[r*cos(a) r*sin(a) xy(i,3)];
end
fens.xyz=xy;

femm = femm_deformation_nonlinear(struct ('material',mater, 'fes',fes, ...
    'integration_rule',gauss_rule(struct('dim',3,'order',2))));

% Material
% Finite element block

efemm = femm_deformation (struct ('material',mater, 'fes',subset(bdry_fes,icl),...
    'integration_rule',surface_integration_rule));
sfemm = femm_deformation (struct ('material',mater, 'fes',bdry_fes,...
    'integration_rule',surface_integration_rule));

% Geometry
geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
% Define the displacement field
u   = 0*geom; % zero out
% Apply EBC's
                node_list = fenode_select(fens,struct('box',[-100*radius 100*radius 0 0 -100*radius 100*radius],'inflate',htol));
u= set_ebc(u,node_list,true,1,0.0);
u= set_ebc(u,node_list,true,2,0.0);
u= set_ebc(u,node_list,true,3,0.0);
                u   = apply_ebc (u);
gv=drawmesh({fens,fes},'fes');
show_field_as_markers(gv,struct('x',geom,'u',u,'nl',node_list,'color','r','markersize',8))
show_field_as_markers(gv,struct('x',geom,'u',u,'nl',connected_nodes(efemm.fes),'color','b','markersize',8))

% Number equations
u   = numberdofs (u);
% Now comes the nonlinear solution
tup = 1;
u = u*0; % zero out the displacement
utol =         utol*u.nfreedofs;
us={};

if (graphics),
    gv=reset(clear(graphic_viewer,[]),[]);
    cmap = jet;
    Cam= [56.1782   -9.9204   62.2090   60.7439   -3.9703   57.8788         0         0    1.0000  151.2154];
end

% Update the FEMM
femm  =associate_geometry(femm,geom); 
incr=1;
while (incr <= nincr)
   lambda = incr* tup / nincr;
    disp(['Increment ' num2str(incr) ]); % pause
    % Initialization
    u1 = u; % guess
    u1 = apply_ebc(u1);
    du = 0*u; % this will hold displacement increment
    du = apply_ebc(du);
    iter=1;
    %             gv=reset(clear(gv,[]),[]);
    %             draw(sfemm,gv, struct ('x', geom,'u',u, 'facecolor','red'));
    femm1=femm;
    while 1
        Load=zeros(3,1); Load(3) =Fmag/h/b*lambda;
        fi=force_intensity(struct('magn',Load));
        FL = distrib_loads(efemm, sysvec_assembler, geom, 0*u, fi, 2);
        F = FL + restoring_force(femm1,sysvec_assembler, geom,u1,u);       % Internal forces
        K = stiffness(femm1, sysmat_assembler_sparse, geom, u1,u) + stiffness_geo(femm1, sysmat_assembler_sparse, geom, u1,u);
        % Displacement increment
        du = scatter_sysvec(du, K\F);
        R0 = dot(F,gather_sysvec(du));
        F = FL + restoring_force(femm1,sysvec_assembler, geom,u1+du,u);       % Internal forces
        R1 = dot(F,gather_sysvec(du));
        a = R0/R1;
        if ( a<0 )
            eta = a/2 +sqrt((a/2)^2 -a);
        else
            eta =a/2;
        end
        if (imag(eta)~=0)
            disp('######################  Inverted elements?')
        end
        eta=min( [eta, 1.0] );
        u1 = u1 + eta*du;   % increment displacement
        disp(['   It. ' num2str(iter) ': ||du||=' num2str(norm(du))]);
        if (max(abs(du.values)) < utol) break; end;                    % convergence check
        iter=iter+1;
    end
    % Update the FEMM
   [~,femm1]  =restoring_force(femm1,sysvec_assembler,geom,u1,u);
   disp(['    Converged for lambda=' num2str(lambda)]); % pause
    u = u1;                                               % update the displacement
    if graphics
        gv=reset(clear(gv,[]),[]);
        draw(sfemm,gv, struct ('x', geom, 'u', 0*u,'facecolor','none', 'shrink',1.0));
        draw(sfemm,gv, struct ('x', geom, 'u', u,'facecolor','y', 'shrink',1.0));
        camset (gv,Cam);
        interact(gv);
        pause(0.5); Cam =camget(gv);
    end
    us{end+1} =u;
    uend=mean(gather_values(u,connected_nodes(efemm.fes)))
    incr = incr + 1;
end

end
