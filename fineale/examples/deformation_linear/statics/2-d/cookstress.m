% Cook membrane problem, plane stress.%
function [h,utip]=cookstress
t=timetic; tic(t);
    E=1;
    nu=1/3;
    width =48; height = 44; thickness  = 1.0;
    free_height  = 16;
    Mid_edge  = [48, 52];% Location of tracked  deflection
    magn=1/free_height;% Magnitude of applied load
    convutip=23.97;
    n=round(sqrt(170)); % number of elements per side
    tolerance=min(width,height)/n/1000;%Geometrical tolerance
    smult=0;% multiplier of random shift of nodes to produce poorly shaped elements
    graphics=true;
    scale=1;
    
    [fens,fes] = T3_block (width,height, n, n, thickness);
    dxy=min(width,height)/n/10;
    sxy=smult*dxy;
    rand('state',[0.3085,0.4953,0.0143,0.3137,0.7750,0.8827,0.6275,0.5996,0.3557,0.8033,0.4425,0.3749,0.3086,0.6245,0.0244,0.0309,0.1962,0.2670,0.8672,0.8259,0.3590,0.6446,0.3018,0.6694,0.5783,0.3251,0.0024,0.9082,0.4464,0.0331,0.9344,0.0261,0,0.0000,0.0000]');
    xy=fens.xyz;
    for i=1:count(fens)
        %             Random shift of internal nodes
        if (xy(i,1)>tolerance) && (xy(i,2)>tolerance) && (xy(i,1)<width-tolerance) && (xy(i,2)<height-tolerance)
            xy(i,:)=[xy(i,1)+(2*rand(1)-1)/2*sxy,xy(i,2)+(2*rand(1)-1)/2*sxy];
        end
        xy(i,:)=[xy(i,1),xy(i,2) + (xy(i,1)/width)*(height -xy(i,2)/height*(height-free_height))];
    end
    fens.xyz =xy;
    
    % Compose the model data
    clear model_data
    model_data.fens =fens;
    
    clear region
    region.E =E;
    region.nu =nu;
    region.reduction ='stress';
    region.fes= fes;
    region.integration_rule = tri_rule (struct('npts', 1));
    region.Rm =[];
    model_data.region{1} =region;
    
    clear essential
    essential.component= [];% clamped
    essential.fixed_value= 0;
    essential.node_list = fenode_select (fens,struct('box',[0,0,-inf, inf],'inflate',tolerance));
    model_data.boundary_conditions.essential{1} = essential;
    
    boundary_fes =  mesh_boundary (fes, struct( 'other_dimension',  thickness ));
    Toplist  =fe_select(fens,boundary_fes,struct('box',  [width, width, -inf, inf ], 'inflate',  tolerance));
    clear traction
    traction.fes= subset(boundary_fes,Toplist);
    traction.integration_rule = trapezoidal_rule (struct('dim', 1));
    traction.traction = [0;magn];
    model_data.boundary_conditions.traction{1} = traction;
    
    
    % Solve
    model_data =deformation_linear_statics(model_data);
    nl=fenode_select (fens,struct ('box',[Mid_edge(1),Mid_edge(1),Mid_edge(2),Mid_edge(2)],'inflate',tolerance));
    theutip=gather_values(model_data.u,nl);
    disp (['displacement =' num2str(theutip(2)) ' as compared to converged ' num2str(convutip)])
    
    model_data.postprocessing.u_scale= scale;
    model_data.postprocessing.draw_mesh= ~false;
    model_data.postprocessing.boundary_only= false;
    model_data.postprocessing.cmap= parula;
    model_data=deformation_plot_deformation(model_data)
    view (2);;
    
    count(fens)
    toc(t)