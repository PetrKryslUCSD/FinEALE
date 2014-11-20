% Three-dimensional Elasticity Solution for Uniformly Loaded Cross-ply
% Laminates and Sandwich Plates
% Ashraf M. Zenkour, Journal of Sandwich Structures and Materials 2007 9: 213-238
% DOI: 10.1177/1099636207065675
%
function Z_laminate_u_ss
    u= physical_units_struct;
    % E1=172350; E2=6894; E3=6894; G12=3447;  G13=3447; G23=1378.8;
    E1=25*u.MEGA*u.PSI; E2=1*u.MEGA*u.PSI; E3=E2; G12=0.5*u.MEGA*u.PSI;  G13=G12; G23=0.2*u.MEGA*u.PSI;
    nu12= 0.25; nu13= 0.25; nu23= 0.25;
    a=200*u.MM; b=600*u.MM;
    q0 = -1*u.PSI;
    % The below values come from Table 2
    h=a/4; wc_analytical  =3.65511/(100*E3*h^3/a^4/q0);
    % h=a/10; wc_analytical  =1.16899/(100*E3*h^3/a^4/q0);
    h=a/50; wc_analytical  =0.66675/(100*E3*h^3/a^4/q0);
    h=a/100; wc_analytical  =0.65071/(100*E3*h^3/a^4/q0);
    
    angles =[0,90,0];
    nLayers =length(angles);
    [na,nb] =adeal(2*[4,4]);
    nts= 1*ones(nLayers,1);% number of elements per layer
    ts= h/nLayers*ones(nLayers,1);% layer thicknesses
    tolerance =min(ts)/max(nts)/100;
    u_scale  = 100;
    graphics = false;
    
    function Rm = LayerRm(XYZ, ts, label)% label equals the layer number here
        Rm= rotmat(angles(label)/180*pi* [0,0,1]);
    end
    % Mesh
    [fens,fes] = H8_composite_plate(a,b,ts,na,nb,nts);;
    [fens,fes] = H8_to_H20(fens,fes);
    % gv=drawmesh( {fens,subset(fes,fe_select(fens,fes,struct('label', 1)))},'fes', 'facecolor','r');
    % gv=drawmesh( {fens,subset(fes,fe_select(fens,fes,struct('label', 2)))},'gv',gv,'fes', 'facecolor','g');
    % gv=drawmesh( {fens,subset(fes,fe_select(fens,fes,struct('label', 3)))},'gv',gv,'fes', 'facecolor','b');
    % gv=drawmesh( {fens,subset(fes,fe_select(fens,fes,struct('label', 4)))},'gv',gv,'fes', 'facecolor','m');
    
    
    % Compose the model data
    clear model_data
    model_data.fens =fens;
    
    clear region
    region.property = 'orthotropic';
    region.E1 =E1;
    region.E2 =E2;
    region.E3 =E3;
    region.G12=G12;
    region.G13=G13;
    region.G23=G23;
    region.nu12=nu12;
    region.nu13=nu13;
    region.nu23=nu23;
    region.fes= fes;
    region.integration_rule = gauss_rule (struct('dim', 3, 'order', 2));
    region.Rm =@LayerRm;
    model_data.region{1} =region;
    
    clear essential
    essential.component= [2,3];
    essential.fixed_value= 0;
    essential.node_list = [[fenode_select(fens, struct('box', [0,0,-Inf,Inf,-Inf,Inf],...
        'inflate',tolerance))],[fenode_select(fens, struct('box', [a,a,-Inf,Inf,-Inf,Inf],...
        'inflate',tolerance))]];;
    model_data.boundary_conditions.essential{1} = essential;
    
    clear essential
    essential.component= [1,3];
    essential.fixed_value= 0;
    essential.node_list = [[fenode_select(fens, struct('box', [0,a,0,0,-Inf,Inf],...
        'inflate',tolerance))],[fenode_select(fens, struct('box', [0,a,b,b,-Inf,Inf],...
        'inflate',tolerance))]];;
    model_data.boundary_conditions.essential{2} = essential;
    
    clear traction
    bdry_fes = mesh_boundary(fes, []);
    bcl = fe_select(fens, bdry_fes, ...
        struct ('box',[-Inf,Inf,-Inf,Inf,h,h],'inflate',tolerance));
    traction.fes =subset(bdry_fes,bcl);;
    traction.traction= [0;0; q0];
    traction.integration_rule =gauss_rule (struct('dim', 2, 'order', 2));
    model_data.boundary_conditions.traction{1} = traction;
    
    % Solve
    model_data =deformation_linear_statics(model_data);
     
    model_data.postprocessing.u_scale= u_scale;
    model_data=deformation_plot_deformation(model_data);
    
    nc=[fenode_select(fens, struct('box', [a/2,a/2,b/2,b/2,-Inf,Inf],...
        'inflate',tolerance))];
    Center_displacement=gather_values(model_data.u,nc);
    wc =mean(Center_displacement (:,3) );
    disp(['Center deflection wc/wc_analytical=' num2str(wc/wc_analytical)])
    
    % clear options
    % options.u_scale= scale;
    % options.stress_component= 1;
    % options=deformation_plot_stress(model_data, options)
end
