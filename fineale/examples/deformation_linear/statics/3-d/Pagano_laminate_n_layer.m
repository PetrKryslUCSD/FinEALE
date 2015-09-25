% This example provides three-dimensional finite element model for  the
% transverse shear stress calculations. The problem consists of a one-, two- or
% three-layer plate subjected to a sinusoidal distributed load, as
% described by Pagano (1969). The resulting transverse shear and axial
% stresses through the thickness of the plate are compared to two existing
% analytical solutions by Pagano (1969). The ?rst solution is derived from
% classical laminated plate theory (CPT), while the second is an exact
% solution from linear elasticity theory.
%
% Span to the thickness ratio changes the effect of shear deformation.
%     One layer [0]
%     For large aspect ratio, the analytical and the CPT solutions give
%     w= 0.494, for aspect ratio=4.0 the elasticity solution yields w=1.969
%     Two layers [0/90]
%     For large aspect ratio, the analytical and the CPT solutions give
%     w= 2.622, for aspect ratio=4.0 the elasticity solution yields w=4.694
%     Three layers [0/90/0]
%     For large aspect ratio, the analytical and the CPT solutions give
%     w= 0.514, for aspect ratio=4.0 the elasticity solution yields w=2.906
    
function Pagano_laminate_n_layer
    u= physical_units_struct;
    E1=25e6*u.PSI; E2=1e6*u.PSI; E3=E2; G12=0.5e6*u.PSI;  G13=G12; G23=0.2e6*u.PSI;
    nu12= 0.25; nu13= 0.25; nu23= 0.25;
    Span_to_thickness=4.0;
   T=2.5*u.IN; L=Span_to_thickness*T; h=1*u.IN; 
    q0 = 1*u.PSI;
    
    %   angles =[0,90,0];% Three-layer plate
          angles =[0,90];% d%bTwo-layer plate
        %         angles =[0];% Single layer plate
    nLayers =length(angles);
    [nL,nh] =adeal([28,1]);
    nts= 8*ones(nLayers,1);% number of elements per layer
    ts= T/nLayers*ones(nLayers,1);% layer thicknesses
    tolerance =min(ts)/max(nts)/100;
    graphics = false;
    
    function Rm = LayerRm(XYZ, ts, label)% label equals the layer number here
        Rm= rotmat(angles(label)/180*pi* [0,0,1]);
    end
    % Mesh
    [fens,fes] = H8_composite_plate(L,h,ts,nL,nh,nts);;
    [fens,fes] = H8_to_H20(fens,fes);
    %     gv=drawmesh( {fens,subset(fes,fe_select(fens,fes,struct('label', 1)))},'fes', 'facecolor','r');
    %     gv=drawmesh( {fens,subset(fes,fe_select(fens,fes,struct('label', 2)))},'gv',gv,'fes', 'facecolor','g');
    % %     % gv=drawmesh( {fens,subset(fes,fe_select(fens,fes,struct('label', 3)))},'gv',gv,'fes', 'facecolor','b');
    % %     % gv=drawmesh( {fens,subset(fes,fe_select(fens,fes,struct('label', 4)))},'gv',gv,'fes', 'facecolor','m');
    %     labels
    
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
    
    %     The bending occurs in the plane  X-Z, out of this plane the deformation is zero
    clear essential
    essential.component= [2];
    essential.fixed_value= 0;
    essential.node_list = (1:count(fens));
    model_data.boundary_conditions.essential{1} = essential;
    
    %     Simply-supported ends
    clear essential
    essential.component= [3];
    essential.fixed_value= 0;
    essential.node_list = [...
        [fenode_select(fens, struct('box', [0,0,0,h,-Inf,Inf],...
        'inflate',tolerance))],...
        [fenode_select(fens, struct('box', [L,L,0,h,-Inf,Inf],...
        'inflate',tolerance))]];;
    model_data.boundary_conditions.essential{2} = essential;
    
    %     Fixed in the axial direction (These points are located on the symmetry plane)
    clear essential
    essential.component= [1];
    essential.fixed_value= 0;
    essential.node_list = [...
        [fenode_select(fens, struct('box', [L/2,L/2,0,h,-Inf,Inf],...
        'inflate',tolerance))]];;
    model_data.boundary_conditions.essential{3} = essential;
    
    clear traction
    bdry_fes = mesh_boundary(fes, []);
    bcl = fe_select(fens, bdry_fes, ...
        struct ('box',[-Inf,Inf,-Inf,Inf,T,T],'inflate',tolerance));
    traction.fes =subset(bdry_fes,bcl);;
    traction.traction= @(x)[0;0; -q0*sin(pi*x(1)/L)];
    traction.integration_rule =gauss_rule (struct('dim', 2, 'order', 4));
    model_data.boundary_conditions.traction{1} = traction;
    
    % Solve
    model_data =deformation_linear_statics(model_data);
     
    nc=[fenode_select(fens, struct('box', [L/2,L/2,0,h,0,T],...
        'inflate',tolerance))];
    Center_displacement=gather_values(model_data.u,nc);
    wc =mean(Center_displacement (:,3) );
    disp(['center vertical displacement =',num2str(wc/u.IN),'[in]'])
    disp(['w parameter =' num2str(abs((100*E2*T^3)/(q0*L^4)*wc))])
    
    model_data.postprocessing.u_scale= abs(L*((100*E2*T^3)/(q0*L^4)))/50;
    model_data.postprocessing.show_mesh= true;
    model_data.postprocessing.quantity= 'U3';
   model_data=deformation_plot_deformation(model_data);
    
    pause(2)
    clear( model_data.postprocessing.gv);
    reset( model_data.postprocessing.gv,[]);
    model_data.postprocessing.u_scale= 0;
    model_data.postprocessing.stress_component= 5;
    model_data=deformation_plot_stress(model_data);
end
