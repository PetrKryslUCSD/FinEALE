% Nonlinear twisted beam.
%

function [u1s,lambdas]=nonlinear_twisted_beam
E=0.29e8;
nu=0.22;
W=1.1;
L=12;
t= 0.05;
nl=4; nt=2; nw=2;
ref=2;
p=  1/W/t;
loadv=[0;0;p];dir=3;uzex=0.005424534868469; % Harder: 5.424e-3;
loadv=[0;p;0];dir=2;uzex=0.001753248285256; % Harder: 1.754e-3;
utol = t/1e10;
graphics = true;
nincr  =10;

% prop = property_deformation_neohookean (struct('E',E,'nu',nu));
% mater = material_deformation_neohookean_triax(struct('property',prop));
prop = property_deformation_linear_iso (struct('E',E,'nu',nu));
mater = material_deformation_stvk_triax(struct('property',prop));

surface_integration_rule=gauss_rule(struct('dim',2,'order',2));


%          Mesh
%% Create the mesh and initialize the geometry
[fens,fes]= H8_block(L,W,t, nl*ref,nw*ref,nt*ref);
xy=fens.xyz;
for i=1:count (fens)
    a=xy(i,1)/L*(pi/2); y=xy(i,2)-(W/2); z=xy(i,3)-(t/2);
    xy(i,:)=[xy(i,1),y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
end
fens.xyz=xy;

femm = femm_deformation_nonlinear_h8msgso(struct ('material',mater, 'fes',fes, ...
    'integration_rule',gauss_rule(struct('dim',3,'order',2))));

% Material
% Finite element block
bdry_fes = mesh_boundary(fes, []);
            bcl = fe_select(fens, bdry_fes, ...
                struct ('box',[L L -100*W 100*W -100*W 100*W],'inflate',0.0001*t));
            
efemm = femm_deformation (struct ('material',mater, 'fes',subset(bdry_fes,bcl),...
    'integration_rule',surface_integration_rule));
sfemm = femm_deformation (struct ('material',mater, 'fes',bdry_fes,...
    'integration_rule',surface_integration_rule));
% Geometry
geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
% Define the displacement field
u   = 0*geom; % zero out
% Apply EBC's
fixed_value= 0;
node_list = fenode_select (fens,struct ('box',[0 0 -100*W 100*W -100*W 100*W],'inflate',0.001*t));;
u= set_ebc(u,node_list,true,1,fixed_value);
u= set_ebc(u,node_list,true,2,fixed_value);
u= set_ebc(u,node_list,true,3,fixed_value);
u   = apply_ebc (u);

% Number equations
u   = numberdofs (u);
% Now comes the nonlinear solution
tup = 60;
u = u*0; % zero out the displacement
utol =         utol*u.nfreedofs;
us={};

if (graphics),
    gv=reset(clear(graphic_viewer,[]),[]);
    cmap = jet;
    Cam= [-0.171332705794298  -7.882139089645855   5.594516542362686   4.394378011771107  -1.931989037461593   1.264389523440495                   0   0   1.000000000000000  54.988185976473318];
end

enl=fenode_select (fens,struct ('box',[L L -100*W 100*W -100*W 100*W],'inflate',0.01*t));
            

u1s=[]; lambdas=[];
% Update the FEMM
femm  =associate_geometry(femm,geom); 
incr=1;
while (incr <= nincr)
    t = incr* tup / nincr;
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
        Load=zeros(3,1); Load(2) =p*t;
        fi=force_intensity(struct('magn',Load));
        FL = distrib_loads(efemm, sysvec_assembler, geom, 0*u, fi, 2);
        %         sum(FL)
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
    [~,femm]  =restoring_force(femm1,sysvec_assembler,geom,u1,u);
    disp(['    Converged for t=' num2str(t)]); % pause
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
    uend=mean(gather_values(u,enl));
    u1s=[u1s,uend];
    lambdas=[lambdas,t];
    incr = incr + 1;
end
u1s  =reshape(u1s',[3,nincr])';
end
