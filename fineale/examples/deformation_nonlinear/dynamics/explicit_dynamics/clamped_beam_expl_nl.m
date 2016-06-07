function clamped_beam_expl_nl
    % This example concerns the response of an elastic beam, built-in at
    % both ends, subject to a suddenly applied load at its midspan
    %     The central part of the beam undergoes displacements
    % several times its thickness, so the solution quickly becomes
    % dominated by membrane effects that significantly stiffen its
    % response.
    %
    % The double cantilever beam has a span of 508 mm (20 in), with a
    % rectangular cross-section 25.4 mm (1 in) wide by 3.175 mm (0.125 in)
    % deep. The material is linear elastic, with a Young's modulus of 206.8
    % GPa (30 x 10^6 lb/in^2) and a density of 2710.42 kg/m3 (2.5362 x
    % 10^-4 lb-s^2/in^4).
    %
    %     The benchmark "Double cantilever elastic beam under point load" from the Abaqus Benchmarks Guide.
    
    graphics = false;
    igraphics=(1:100:40000000);
    plots = true;
    % Parameters:
    u=physical_units_struct;
    E=206800*u.MEGA*u.PA;
    nu=0.3;
    rho  = 2714.5*u.KG/u.M^3;
    L = 508/2*u.MM;
    H = 3.175*u.MM;
    W = 25.4/2*u.MM;
    Force= 2846.7/4*u.NT;% We are modeling just one quarter
    tolerance  =H/1000;
    tend = 5e-3*u.SEC;
    scale=4000;
    
    prop = property_deformation_linear_iso (struct('E',E,'nu',nu,'rho',rho));
    mater = material_deformation_stvk_triax(struct('property',prop));
    
    %     Mesh
    if (1)% Choose hexahedral mesh
        [fens,fes] = H8_block(L,W,H, 4,1,2);
        femm = femm_deformation_nonlinear_h8msgso(struct ('material',mater, 'fes',fes, ...
            'integration_rule',gauss_rule(struct('dim',3,'order',2))));
        Surface_integration_rule =gauss_rule(struct('dim',2, 'order',2));
    else% Choose tetrahedral mesh
        [fens,fes] = T10_blockb(L,W,H, 4,1,2);
        [fens,fes] = T10_to_T10MS(fens,fes);
        femm = femm_deformation_nonlinear_t10ms(struct ('material',mater, 'fes',fes, ...
            'integration_rule',tet_rule(struct('npts',1))));
        Surface_integration_rule =tri_rule(struct('npts',1));
    end
    
    
    % Compose the model data
    clear model_data
    model_data.fens =fens;
    
    
    
    clear region
    region.femm = femm;
    model_data.region{1} =region;
    
    
    clear essential
    essential.component= [1];% symmetry with respect to the YZ plane
    essential.fixed_value= 0;
    essential.node_list = fenode_select (fens,struct('box',[0,0,-Inf,Inf,-Inf,Inf],'inflate',tolerance));
    model_data.boundary_conditions.essential{1} = essential;
    clear essential
    essential.component= [2];% symmetry with respect to the XZ plane
    essential.fixed_value= 0;
    essential.node_list = fenode_select (fens,struct('box',[-Inf,Inf,0,0,-Inf,Inf],'inflate',tolerance));
    model_data.boundary_conditions.essential{2} = essential;
    clear essential
    essential.component= [];% clamped end
    essential.fixed_value= 0;
    essential.node_list = fenode_select (fens,struct('box',[L,L,-Inf,Inf,-Inf,Inf],'inflate',tolerance));
    model_data.boundary_conditions.essential{3} = essential;
    
    %%
    % Define the traction load of surface of the symmetry cut.  
    clear traction
    traction.traction = [0;0;Force/W/H];
    bfes= mesh_boundary(fes,[]);
    topl =fe_select (fens,bfes,struct('facing', true,'direction', [-1,0,0]));
    traction.fes= subset(bfes,topl);
    traction.integration_rule =Surface_integration_rule;
    model_data.boundary_conditions.traction{1} = traction;
    
    
    clear initial_condition
    initial_condition.u_fixed_value= @(x)0*x;
    initial_condition.v_fixed_value= @(x)0*x;
    model_data.initial_condition = initial_condition;
    
    model_data.tend = tend;;
    model_data.observer  =@output;;
    midp=fenode_select (fens,struct('box',[0 0 0 0 0 0],'inflate',tolerance));
    midpu  = [];
    if graphics
        gv=graphic_viewer;
        gv=reset (gv,[]);
        peekb  = uicontrol('Parent',gcf,'Style','pushbutton',...
                    'Units','Normalized','Position',[.9 .0 .1 .05],...
                    'String','Peek','Value',0,...
                    'TooltipString','Invoke command line in order to peek at data',...
                    'Callback',@(h,varargin)keyboard);
    elseif plots
        pf = figure(gcf);
    end
    
    % Solve
%         model_data.dt=1e-8;
    tic; model_data =deformation_nonlinear_direct_explicit_CD(model_data);
    toc
    
    function output(t, model_data)
        midpu= [midpu    model_data.un1.reshape(gather_values(model_data.un1, midp))];
        if (~isempty(igraphics) && (igraphics(1)==round(t/model_data.dt)))
            figure (pf);
            plot(t/(u.MILLI*u.SEC),midpu(3,end)/u.IN,'rx','Markersize',4); hold on
            labels ('Time [ms]','Deflection [in]')
            pause(0.01)
            igraphics= igraphics(2:end);
        end
        %         end
    end
end