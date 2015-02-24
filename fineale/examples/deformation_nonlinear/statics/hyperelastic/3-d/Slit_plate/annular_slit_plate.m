% Annular shell (slit plate)
%
function [lambdas,uAs,uBs]=annular_slit_plate
u= physical_units_machine;
E=21e6;%*u('Pa')
nu=0.;
R= 6;%*u('mm');
Width=4;%*u('mm');;%
Thickness=0.03;%*u('mm');
ang=360/180*pi;
p=  0.8/Thickness;
nc=80;nT=4;nW=8;
nc=40;nT=3;nW=4;
nc=30;nT=3;nW=6;
nc=40;nT=3;nW=8;
nc=20;nT=2;nW=8;
nc=20;nT=2;nW=3;

nincr = 10;
graphics = ~false;
scale=1;
utol = Thickness/1e6;

prop = property_deformation_neohookean (struct('E',E,'nu',nu));
mater = material_deformation_neohookean_triax(struct('property',prop));

surface_integration_rule=gauss_rule(struct('dim',2,'order',2));


%          Mesh
%% Create the mesh and initialize the geometry
[fens,fes]= H8_block(ang, Thickness, Width,nc,nT,nW);
bdry_fes = mesh_boundary(fes, struct('other_dimension',1.0));
l1=fe_select(fens,bdry_fes,struct('box',[0 0 -inf inf -inf inf],'inflate',Thickness/100));
l2=fe_select(fens,bdry_fes,struct('box',[ang ang -inf inf -inf inf],'inflate',Thickness/100));
Al=fenode_select(fens,struct('box',[ang ang -inf inf 0 0],'inflate',Thickness/100));
Al=intersect(connected_nodes(subset(bdry_fes,l2)),Al);
Bl=fenode_select(fens,struct('box',[ang ang -inf inf Width Width],'inflate',Thickness/100));
Bl=intersect(connected_nodes(subset(bdry_fes,l2)),Bl);

xy=fens.xyz;
for i=1:count (fens)
    a=xy(i,1); y=xy(i,2); r=R+xy(i,3);
    xy(i,:)=[r*sin(a) y (r*cos(a))];
end
fens.xyz=xy;

femm = femm_deformation_nonlinear_h8msgso(struct ('material',mater, 'fes',fes, ...
    'integration_rule',gauss_rule(struct('dim',3,'order',2))));

% Material
% Finite element block
efemm = femm_deformation (struct ('material',mater, 'fes',subset(bdry_fes,l2),...
    'integration_rule',surface_integration_rule));
sfemm = femm_deformation (struct ('material',mater, 'fes',bdry_fes,...
    'integration_rule',surface_integration_rule));
% Geometry
geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
% Define the displacement field
u   = 0*geom; % zero out
% Apply EBC's
fixed_value= 0;
node_list = connected_nodes (subset(bdry_fes,l1));
u= set_ebc(u,node_list,true,1,fixed_value);
u= set_ebc(u,node_list,true,2,fixed_value);
u= set_ebc(u,node_list,true,3,fixed_value);
u   = apply_ebc (u);

% Number equations
u   = numberdofs (u);
% Now comes the nonlinear solution
tup = 1;
u = u*0; % zero out the displacement
% utol =         utol*u.nfreedofs;
us={};

if (graphics),
    gv=reset(clear(graphic_viewer,[]),[]);
    cmap = jet;
    Cam= 1.0e+02 *[-0.026964096549535  -0.052147729718727   0.061678333777267   0.005000000000000   0.005000000000000   0.005000000000000   0.002021395347527  0.006303304978612   0.007495485787833   1.004597304098127];
    Cam= [ -4.065710717565406  -5.450150052184264   4.830127018922192   0.500000000000000   0.500000000000000   0.500000000000000                   0   0   1.000000000000000  90.095728615633746];
end

lambdas = []; uAs = []; uBs = []; niters=[];
femm  =update(femm,geom,u,u);
incr=1;
while (incr <= nincr)
    t = incr* tup / nincr;
    disp(['Increment ' num2str(incr) ]); % pause
    % Initialization
    u1 = u; % guess
    u1 = apply_ebc(u1);
    du = 0*u; % this will hold displacement increment
    du = apply_ebc(du);
    %             gv=reset(clear(gv,[]),[]);
    %             draw(sfemm,gv, struct ('x', geom,'u',u, 'facecolor','red'));
    femm1=femm;
    prevndu=inf;
    iter=1;
    while 1
        Load=zeros(3,1); Load(2) =p*t;
        fi=force_intensity(struct('magn',Load));
        FL = distrib_loads(efemm, sysvec_assembler, geom, 0*u, fi, 2);
        F = FL + restoring_force(femm1,sysvec_assembler, geom,u1,u);       % Internal forces
        K = stiffness(femm1, sysmat_assembler_sparse, geom, u1,u) + stiffness_geo(femm1, sysmat_assembler_sparse, geom, u1,u);
        % Displacement increment
        du = scatter_sysvec(du, K\F);
        ndu=norm(du);;
        if (iter<10)||(ndu>prevndu/2)
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
        else
            eta= 1.0;
        end
        u1 = u1 + eta*du;   % increment displacement
        disp(['   It. ' num2str(iter) ': ||du||=' num2str(ndu) ' (eta=' num2str(eta) ')']);
        if (ndu < utol) break; end;  % convergence check
        prevndu =ndu;
        iter=iter+1;
    end
    femm  =update(femm,geom,u1,u);
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
    uA=mean(gather_values(u,Al));
    uB=mean(gather_values(u,Bl));;
    lambdas = [lambdas,t]; uAs = [uAs,uA]; uBs = [uBs,uB];
    niters=[niters,iter];
    incr = incr + 1;
end
save('res8','lambdas','uAs','uBs','niters');
end
