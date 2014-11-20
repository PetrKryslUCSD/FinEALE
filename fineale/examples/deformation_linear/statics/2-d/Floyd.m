% Floyd's pressure vessel example (from the book by Babuska, Szabo)
function  Floyd
    u=  physical_units_struct;
    % Parameters:
    E= 1435*u.PSI;% 
    nu=0.49;
    % geometry
    R1= 6*u.IN;%  inches
    Rf=0.15*u.IN;% Inches
    R2=R1-0.6*u.IN;% inches
    L1=3.15*u.IN;%Inches
    L2=1.5*u.IN;%Inches
    % prescribed Traction
    Pressure=2.61*u.PSI;%PSI
    sigj=3;
    scale=5;

    % Mesh'
    mesh_size=R1/(2^2);
    [fens,fes,groups,edge_fes,edge_groups]=targe2_mesher({...
        ['curve 1 line ' num2str(0) ' ' num2str(0) ' ' num2str(R1) ' ' num2str(0) ],...
        ['curve 2 line ' num2str(R1) ' ' num2str(0) ' ' num2str(R1) ' ' num2str(L1) ],...
        ['curve 3 line ' num2str(R1) ' ' num2str(L1) ' ' num2str(R2) ' ' num2str(L1) ],...
        ['curve 4 line ' num2str(R2) ' ' num2str(L1) ' ' num2str(R2) ' ' num2str(L2+Rf) ],...
        ['curve 5 arc ' num2str(R2) ' ' num2str(L2+Rf) ' ' num2str(R2-Rf) ' ' num2str(L2) ' center ' num2str(R2-Rf) ' ' num2str(L2+Rf) ],...
        ['curve 6 line ' num2str(R2-Rf) ' ' num2str(L2) ' ' num2str(0) ' ' num2str(L2) ],...
        ['curve 7 line ' num2str(0) ' ' num2str(L2) ' ' num2str(0) ' ' num2str(0) ],...
        ['subregion 1  property 1 boundary 1 2 3 4 5 6 7'],...
        ['m-ctl-point constant ' num2str(mesh_size)],...
        ['m-ctl-point 1 xy ' num2str(R2-Rf) ' ' num2str(L2+Rf) ' near ' num2str(mesh_size/20000) ' influence ' num2str(Rf/100)],...
        }, 1.0, struct('axisymm',true,'quadratic',true));
    %     drawmesh({fens,fes},'fes','facecolor','red');  view(2); return


    % Compose the model data
    clear model_data
    model_data.fens =fens;
    
    clear region
    region.E =E;
    region.nu =nu;
    region.reduction ='axisymm';
    region.fes= fes;
    region.integration_rule = tri_rule (struct('npts', 3));
    region.Rm =[];
    model_data.region{1} =region;
    
    clear essential % the symmetry plane
    essential.component= 2;
    essential.fixed_value= 0;
    essential.node_list = fenode_select (fens,struct('box',[0 R1 L1 L1],'inflate',L1/1000));
    model_data.boundary_conditions.essential{1} = essential;
    
    clear essential % the axis of symmetry
    essential.component= 1;
    essential.fixed_value= 0;
    essential.node_list = fenode_select (fens,struct('box',[0 0 0 L1+L2],'inflate',L1/1000));
    model_data.boundary_conditions.essential{2} = essential;
    
    clear traction
    traction.fes= subset(edge_fes,edge_groups{4});
    traction.integration_rule = gauss_rule (struct('dim', 1,  'order', 2));
    traction.traction = [Pressure;0];
    model_data.boundary_conditions.traction{1} = traction;
    
    clear traction
    traction.fes= subset(edge_fes,edge_groups{5});
    traction.integration_rule = gauss_rule (struct('dim', 1,  'order', 2));
    traction.traction = @(x) (Pressure*(x-[R2-Rf,L2+Rf])'/norm((x-[R2-Rf,L2+Rf])));
    model_data.boundary_conditions.traction{2} = traction;
    
    clear traction
    traction.fes= subset(edge_fes,edge_groups{6});
    traction.integration_rule = gauss_rule (struct('dim', 1,  'order', 2));
    traction.traction = [0;-Pressure];
    model_data.boundary_conditions.traction{3} = traction;
    
    
    % Solve
    model_data =deformation_linear_statics(model_data);
    
    model_data.postprocessing.u_scale= scale;
    model_data.postprocessing.stress_component=sigj;
    model_data.postprocessing.stress_range = 10*[-Pressure,+Pressure];
    model_data=deformation_plot_stress(model_data);
    view (2);;
    
    