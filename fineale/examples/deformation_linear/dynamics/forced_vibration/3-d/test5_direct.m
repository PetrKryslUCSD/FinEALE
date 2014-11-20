% Construction of the frequency response function by direct integration.
% This is a test recommended by the National Agency for Finite Element
% Methods and Standards (U.K.): Test 5 from NAFEMS “Selected Benchmarks
% for Forced Vibration,” R0016, March 1993.
%
% The extensional and torsional modes are not correct due to the
% impossibility of enforcing beam-like point supports, but the flexural
% modes are within a percent of the analytical value.
function test5_direct
    pu=physical_units_struct;
    % Parameters:
    E = 200e3*pu.MEGA*pu.PA;
    nu = 0.3;
    rho= 8000*pu.KG/pu.M^3;
    a=2.00*pu.M; b=2.00*pu.M; L= 10*pu.M;
    tolerance =a/1000;
    frequency = 42.7;% forcing frequency
    zeta= 0.02;% damping ratio
    tend= 60/frequency;
    dt= (1/frequency)/20;
    omega =2*pi*frequency;
    %     eta = loss_tangent*G/(2*pi*frequency);
    %  Reference:   RHEOLOGICAL INTERPRETATION OF RAYLEIGH DAMPING, Semblat
    %         Rayleigh_mass = 2*zeta*(2*pi*frequency);
    %         Rayleigh_stiffness = 2*zeta/(2*pi*frequency);
    zeta1= zeta; zeta2  =zeta;
    f1= 42.75; f2 =  149.4; o1 =2*pi*f1;  o2 =2*pi*f2;
    Rayleigh_mass = 2*(o1*o2)/(o2^2-o1^2)*(o2*zeta1-o1*zeta2);
    Rayleigh_stiffness = 2*(o1*o2)/(o2^2-o1^2)*(-1/o2*zeta1+1/o1*zeta2);
    graphics = ~true;
    igraphics=(1:1:40000000);
    plots = ~graphics;
    scale = 5000;
    % Mesh refinement
    na= 2; nb=  2; nL =10;
    %     na= 4; nb=  4; nL =20;
    
    prop=property_deformation_linear_iso(struct('E',E,'nu',nu,'rho',rho));
    mater = material_deformation_linear_triax (struct('property',prop ));
    
    clear eltyd
    eix=1;
    
    %         eltyd(eix).description ='H64';
    %         eltyd(eix).mf =@H64_block;
    %         eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
    %             'integration_rule',gauss_rule(struct('dim',3, 'order',4))));
    %         eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',4));
    %         eltyd(eix).styl='k*--';
    %         eix=eix+1;
    %
    %         eltyd(eix).description ='H27';
    %         eltyd(eix).mf =@H27_block;
    %         eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
    %             'integration_rule',gauss_rule(struct('dim',3, 'order',3))));
    %         eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',4));
    %         eltyd(eix).styl='ks--';
    %         eix=eix+1;
    %
%     eltyd(eix).description ='T10';% tetrahedron
%     eltyd(eix).mf =@T10_block;
%     eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes,'material',mater,'integration_rule',tet_rule(struct('npts',4))));
%     eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
%     eltyd(eix).styl='k^--';
%     eix=eix+1;
    %
    %         eltyd(eix).description ='T10-SRI';
    %         eltyd(eix).mf =@T10_block;
    %         eltyd(eix).femmf =@(fes)femm_deformation_linear_sri(struct('fes',fes, 'material',mater,...
    %             'integration_rule_volumetric',tet_rule(struct('npts',1)),...
    %             'integration_rule_deviatoric',tet_rule(struct('npts',4)),...
    %             'integration_rule',tet_rule(struct('npts',4))));
    %         eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
    %         eltyd(eix).styl='kv-';
    %         eix=eix+1;
    
    %
        eltyd(eix).description ='H20R';
        eltyd(eix).mf =@H20_block;
        eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
            'integration_rule',gauss_rule(struct('dim',3, 'order',3))));
        eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',4));
        eltyd(eix).styl='ro--';
        eix=eix+1;
    %
    %     %         % Selective reduced integration hexahedron
    %     eltyd(eix).description ='H8-SRI';
    %     eltyd(eix).mf =@H8_block;
    %     eltyd(eix).femmf =@(fes)femm_deformation_linear_sri(struct('fes',fes, 'material',mater,...
    %         'integration_rule_volumetric',gauss_rule(struct('dim',3, 'order',1)),...
    %         'integration_rule_deviatoric',gauss_rule(struct('dim',3, 'order',2))));
    %     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    %     eltyd(eix).styl='md--';
    %     eix=eix+1;
    %
    
    
    ns=[];
    uzs=[];
    
    eix=1;
    %% Create the mesh and initialize the geometry
    [fens,fes]= eltyd(eix).mf(L,a,b,nL,na,nb);
    fens = transform_apply(fens,@(x,d)(x-[0,a/2,b/2]), []);
    bfes= mesh_boundary(fes,[]);
    topl =fe_select (fens,bfes,struct('box', [0,L,-Inf,Inf,b/2,b/2],...
        'inflate',tolerance));
    
    %                     drawmesh({fens,fes},'shrink', 0.8,'facecolor','red');
    
    % Compose the model data
    clear model_data
    model_data.fens =fens;
    
    clear region
    region.fes= fes;
    region.femm= eltyd(eix).femmf(fes);
    model_data.region{1} =region;
    
    clear essential
    essential.component= [1];
    essential.fixed_value= 0;
    essential.node_list = [fenode_select(fens, struct('box', [0,0,0,0,0,0],...
        'inflate',tolerance))];
    model_data.boundary_conditions.essential{1} = essential;
    
    clear essential
    essential.component= [2];
    essential.fixed_value= 0;
    essential.node_list = [fenode_select(fens, struct('box', [0,0,-Inf,Inf,-Inf,Inf],...
        'inflate',tolerance))];
    model_data.boundary_conditions.essential{2} = essential;
    
    clear essential
    essential.component= [3];
    essential.fixed_value= 0;
    essential.node_list = [fenode_select(fens, struct('box', [0,0,-Inf,Inf,-Inf,Inf],...
        'inflate',tolerance))];
    model_data.boundary_conditions.essential{3} = essential;
    
    clear essential
    essential.component= [2];
    essential.fixed_value= 0;
    essential.node_list = [fenode_select(fens, struct('box', [L,L,-Inf,Inf,-Inf,Inf],...
        'inflate',tolerance))];
    model_data.boundary_conditions.essential{4} = essential;
    
    clear essential
    essential.component= [3];
    essential.fixed_value= 0;
    essential.node_list = [fenode_select(fens, struct('box', [L,L,-Inf,Inf,-Inf,Inf],...
        'inflate',tolerance))];
    model_data.boundary_conditions.essential{5} = essential;
    
    
    clear traction
    traction.traction = @(t)sin(omega*t)*[0;0;1*pu.MEGA*pu.NT/b];%(1-exp(-t/(2/frequency)))*
    traction.fes= subset(bfes,topl);
    traction.integration_rule =eltyd(eix).surface_integration_rule;
    model_data.boundary_conditions.traction{1} = traction;
    
    clear initial_condition
    initial_condition.u_fixed_value= @(x)0*x;
    initial_condition.v_fixed_value= @(x)0*x+ones(size(x,1),1)*[0,0,0];
    model_data.initial_condition = initial_condition;
    
    model_data.Rayleigh_mass =Rayleigh_mass;
     model_data.Rayleigh_stiffness =Rayleigh_stiffness;
    model_data.tend = tend;
    model_data.dt = dt;
    model_data.observer  =@output;
    midpoint=fenode_select (fens,struct('box',[L/2 L/2 0 0 0 0],'inflate',tolerance));
    midpointu  = [];
    if graphics
        gv=graphic_viewer;
        gv=reset (gv,[]);
    elseif plots
        pf = figure(gcf);
    end
    
    % Solve
    model_data =deformation_linear_direct_implicit_TRAP_Rayleigh(model_data);
    
    function output(t, model_data)
        if graphics && (~isempty(igraphics) && (igraphics(1)==round(t/model_data.dt)))
            %         clear(gv);
            gv=reset (gv,[]);
            draw(model_data.region{1}.fes, gv, struct ('x', model_data.geom, 'u', scale*model_data.u,'facecolor','none','edgecolor','blue'));
            view(3); pause(0.1);
            %             saveas(gcf, [mfilename '-' num2str(etaw) '-' num2str(snapshot) '.png'], 'png');
            igraphics= igraphics(2:end);
        elseif plots
            midpointu= [midpointu    model_data.u.reshape(gather_values(model_data.u, midpoint))];
            if (~isempty(igraphics) && (igraphics(1)==round(t/model_data.dt)))
                figure (pf);
                plot(t,midpointu(3,end),'bx','Markersize',3,'linewidth',2); hold on
                igraphics= igraphics(2:end);
            end
        end
    end
end