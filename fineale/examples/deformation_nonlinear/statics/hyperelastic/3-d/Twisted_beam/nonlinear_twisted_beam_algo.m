% Nonlinear twisted beam
%
% Reference: Krysl, P.,  Mean-strain 8-node Hexahedron with Optimized
% Energy-Sampling Stabilization for Large-strain Deformation, submitted to
% IJNME 18 September 2014.
function [u1s]=nonlinear_twisted_beam_h20_algo
    E=0.29e8;
    nu=0.22;
    W=1.1;
    L=12;
    t= 0.05;
    nl=12; nt=2; nw=2;
    ref=1;
    p=  1/W/t;
    maxdu_tol = t/1e7;
    graphics = true;
    nincr  =10;
    tup = 60;
    
    
    %  Create the mesh and initialize the geometry
    [fens,fes]= H8_block(L,W,t, nl*ref,nw*ref,nt*ref);
    xy=fens.xyz;
    for i=1:count (fens)
        a=xy(i,1)/L*(pi/2); y=xy(i,2)-(W/2); z=xy(i,3)-(t/2);
        xy(i,:)=[xy(i,1),y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
    end
    fens.xyz=xy;
    
    
    
    
    % Package model data
    clear model_data;
    model_data.fens =fens;
    
    
    clear region
    prop = property_deformation_linear_iso (struct('E',E,'nu',nu));
    region.femm= femm_deformation_nonlinear_h8msgso(...
        struct ('material',material_deformation_stvk_triax(struct('property',prop)),...
        'fes',fes, ...
        'integration_rule',gauss_rule(struct('dim',3,'order',2))));;
    model_data.region{1} =region;
    
    %  Clamped cross-section
    clear essential
    essential.component= [1,2,3];
    essential.fixed_value= 0;
    essential.node_list = fenode_select (fens,struct ('box',[0 0 -100*W 100*W -100*W 100*W],'inflate',0.001*t));;
    model_data.boundary_conditions.essential{1} = essential;
    
    % Traction--loaded cross-section
    bdry_fes = mesh_boundary(fes, []);
    bcl = fe_select(fens, bdry_fes, ...
        struct ('box',[L L -100*W 100*W -100*W 100*W],'inflate',0.0001*t));
    clear traction
    traction.fes =subset(bdry_fes,bcl);;
    traction.traction= @(x) ([0,p,0]);
    traction.integration_rule =gauss_rule(struct('dim', 2,'order', 2));
    model_data.boundary_conditions.traction{1} = traction;
    
    % If online graphics  is needed, initialize some variables
    if (graphics),
        sfemm = femm_deformation (struct ('material',[], 'fes',bdry_fes,...
            'integration_rule',[]));
        gv=reset(clear(graphic_viewer,[]),[]);
        cmap = jet;
        Cam= [-0.171332705794298  -7.882139089645855   5.594516542362686   4.394378011771107  -1.931989037461593   1.264389523440495                   0   0   1.000000000000000  54.988185976473318];
    end
    
    % Select the solver options
    model_data.load_multipliers=(1:nincr)/nincr*tup;
    model_data.maxdu_tol  =maxdu_tol;;
    model_data.line_search  = true;
    model_data.iteration_observer =@iteration_observer;
    us={};
    model_data.load_increment_observer =@load_increment_observer;
    % Call the nonlinear deformation solver
    model_data =deformation_nonlinear_statics(model_data);
    
    %     Report results
    enl=fenode_select (fens,struct ('box',[L L -100*W 100*W -100*W 100*W],'inflate',0.01*t));
    u1s=[];
    for j=1:length(us)
        u1s=[u1s,mean(gather_values(us{j},enl))];
    end
    u1s  =reshape(u1s',[3,nincr])';
    
    if graphics
        for j=1:length(us)
            gv=reset(clear(gv,[]),[]);
            draw(sfemm,gv, struct ('x', model_data.geom, 'u', 0*us{j},'facecolor','none', 'shrink',1.0));
            draw(sfemm,gv, struct ('x', model_data.geom, 'u', us{j},'facecolor','y', 'shrink',1.0));
            camset (gv,Cam);
            interact(gv);
            pause(0.05); Cam =camget(gv);
        end
    end
    
    % Observer function to be called when convergence is reached.
    function load_increment_observer(lambda,model_data)
        fprintf(1,'\n');
        if graphics
            gv=reset(clear(gv,[]),[]);
            draw(sfemm,gv, struct ('x', model_data.geom, 'u', 0*model_data.un1,'facecolor','none', 'shrink',1.0));
            draw(sfemm,gv, struct ('x', model_data.geom, 'u', model_data.un1,'facecolor','y', 'shrink',1.0));
            camset (gv,Cam);
            interact(gv);
            pause(0.5); Cam =camget(gv);
        end
        us{end+1} =model_data.un1;
    end
    
    % Iteration of observer can be called as the solution is being computed.
    function iteration_observer(lambda,iter,du,model_data)
        fprintf(1,'+');
        if 0 && graphics
            gv=reset(clear(gv,[]),[]);
            draw(sfemm,gv, struct ('x', model_data.geom, 'u', 0*model_data.un1,'facecolor','none', 'shrink',1.0));
            draw(sfemm,gv, struct ('x', model_data.geom, 'u', model_data.un1,'facecolor','y', 'shrink',1.0));
            camset (gv,Cam);
            interact(gv);
            pause(0.5); Cam =camget(gv);
        end
    end
    
end
