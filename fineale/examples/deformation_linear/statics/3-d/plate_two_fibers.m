% Fiber-reinforced cantilever plate.
function plate_two_fibers
u= physical_units_struct;
% Anisotropic, but much stronger anisotropy
E1=20000*u.MEGA*u.PSI; E2=20000*u.MEGA*u.PSI;  E3=2*u.MEGA*u.PSI; G12=0.5*u.MEGA*u.PSI; G23=0.2*u.MEGA*u.PSI; G13=G23; 
nu12= 0.25; nu13= 0.25; nu23= 0.25;
aangle =45;

% ]; %H8-GSRI
% hs=[2   4   8  16];
% [xestim, beta, c, residual] = richextrapol(uz(end-2:end),hs)
uz_ref =  -5.0733e-006;
prop = property_deformation_linear_ortho (...
    struct('E1',E1,'E2',E2,'E3',E3,...
    'G12',G12,'G13',G13,'G23',G23,...
    'nu12',nu12,'nu13',nu13,'nu23',nu23,...
    'nconstrained',2));


a=90*u.MM; b=100*u.MM;  t=20*u.MM;
q0 = -1*u.PSI;

    refs = [1];
    [na,nb,nt] =adeal(2*[3,3,1]);
tolerance =t/nt/10;
u_scale  = 30000;
graphics = true;

Rm= rotmat(aangle/180*pi* [0,0,1]);

mater = material_deformation_linear_triax (struct('property',prop ));

clear eltyd
eix=1;

% eltyd(eix).description ='THEX';% tetrahedron
% eltyd(eix).mf =@(ref)T15_block_u(a,b,t,2*ref*na,2*ref*nb,2*ref*nt);
% eltyd(eix).femmf =@(fes)femm_deformation_linear_thex(struct('fes',fes,'material',mater,...
%     'Rm',Rm,...
%     'stabfact',0.01*E3));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% eix=eix+1;
%
% eltyd(eix).description ='H20R';
% eltyd(eix).mf =@(ref)H20_block(a,b,t,ref*na,ref*nb,ref*nt);
% eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
%     'Rm',Rm,...
%     'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% eix=eix+1;
% 
%     eltyd(eix).description ='T10';% tetrahedron
%     eltyd(eix).mf =@(ref)T10_block(a,b,t,ref*na,ref*nb,ref*nt);
%     eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes,'material',mater,'Rm',Rm,...
%         'integration_rule',tet_rule(struct('npts',4))));
%     eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
%     eix=eix+1;
%     
% eltyd(eix).description ='NICE-H8';% tetrahedron
% eltyd(eix).mf =@(ref)H8_block(a,b,t,2*ref*na,2*ref*nb,2*ref*nt);
% eltyd(eix).femmf =@(fes)femm_deformation_linear_nice(struct('fes',fes,'material',mater,...
%     'Rm',Rm,...
%     'integration_rule',tensprod_nq_rule(struct('dim',3, 'order',1))));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% eix=eix+1;

% eltyd(eix).description ='H8NIS-GSRI';% tetrahedron
% eltyd(eix).mf =@(ref)H8_block(a,b,t,2*ref*na,2*ref*nb,2*ref*nt);
% eltyd(eix).femmf =@(fes)femm_deformation_linear_nice_1pt_h8_gsri(struct('fes',fes,...
%     'material',mater,...
%     'Rm',Rm));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% eix=eix+1;

%         eltyd(eix).description ='H64';% tetrahedron
%         eltyd(eix).mf =@(ref)H64_block(a,b,t,ref*na,ref*nb,ref*nt);
%         eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes,'material',mater,...
%             'Rm',Rm,...
%             'integration_rule',gauss_rule(struct('dim',3, 'order',4))));
%         eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',3));
%         eix=eix+1;

%                 eltyd(eix).description ='H64-GSRI';
%         eltyd(eix).mf =@(ref)H64_block(a,b,t,ref*na,ref*nb,ref*nt);
%         eltyd(eix).femmf =@(fes)femm_deformation_linear_gsri(struct('fes',fes, 'material',mater,...
%             'Rm',Rm,...
%             'integration_rule_constrained', hex_27pt_rule(),...
%             'integration_rule_unconstrained',hex_64pt_rule()...
%             ));
%         eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',4));
%         eltyd(eix).styl='kx-';
%         eix=eix+1;


% eltyd(eix).description ='H8-GSRI';% tetrahedron
% eltyd(eix).mf =@(ref)H8_block(a,b,t,2*ref*na,2*ref*nb,2*ref*nt);
% eltyd(eix).femmf =@(fes)femm_deformation_linear_gsri(struct('fes',fes,'material',mater,...
%     'Rm',Rm,...
%     'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',1)),...
%     'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',2))));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% eix=eix+1;
% %     %
        eltyd(eix).description ='H27-GSRI';% tetrahedron
        eltyd(eix).mf =@(ref)H27_block(a,b,t,ref*na,ref*nb,ref*nt);
        eltyd(eix).femmf =@(fes)femm_deformation_linear_gsri(struct('fes',fes,'material',mater,...
            'Rm',Rm,...
            'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',2)),...
            'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',3))));
        eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',3));
        eix=eix+1;
        
        % eltyd(eix).description ='H8-Bbar';% tetrahedron
        % eltyd(eix).mf =@(ref)H8_block(a,b,t,2*ref*na,2*ref*nb,2*ref*nt);
        % eltyd(eix).femmf =@(fes)femm_deformation_linear_bbar(struct('fes',fes,'material',mater,...
        %     'Rm',Rm,...
        %     'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',2)),...
        %     'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',2)),...
        %     'pv_bfun',@(p)eye(2)));
        % eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
        % eix=eix+1;
        % %
        %
        % eltyd(eix).description ='H27-Bbar';% tetrahedron
        % eltyd(eix).mf =@(ref)H27_block(a,b,t,ref*na,ref*nb,ref*nt);
        % eltyd(eix).femmf =@(fes)femm_deformation_linear_bbar(struct('fes',fes,'material',mater,...
        %     'Rm',Rm,...
        %     'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',2)),...
        %     'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',3)),...
        %     'pv_bfun',@(p)[eye(2);p(1)*eye(2);p(2)*eye(2);p(3)*eye(2)]));
        % eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',3));
        % eix=eix+1;
%
    %
    % eltyd(eix).description ='H20-GSRI';
    % eltyd(eix).mf =@(ref)H20_block(a,b,t,ref*na,ref*nb,ref*nt);
    % eltyd(eix).femmf =@(fes)femm_deformation_linear_gsri(struct('fes',fes, 'material',mater,...
    %     'Rm',Rm,...
    %     'integration_rule_constrained', hex_8pt_rule(),...
    %     'integration_rule_unconstrained',hex_9pt_rule()...
    %     ));
    % eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',3));

legends ={};
for eix = 1:length(eltyd)
    ns=[];
    uzs=[];
    uys=[];
    ys=[];
    
    disp(['%' eltyd(eix).description ': ' ])
    disp(['uz=[' ])
    %     refs = [2,4,8,16];
    for ref = refs
        
        % Mesh
        [fens,fes] = eltyd(eix).mf(ref);
        
        % Compose the model data
        clear model_data
        model_data.fens =fens;
        
        clear region
        region.fes= fes;
        region.femm= eltyd(eix).femmf(fes);
        model_data.region{1} =region;
        
        clear essential
        essential.component= [1,2,3];
        essential.fixed_value= 0;
        essential.node_list = fenode_select(fens, struct('box', [0,0,-Inf,Inf,-Inf,Inf],...
            'inflate',tolerance));
        model_data.boundary_conditions.essential{1} = essential;
        
        clear traction
        bdry_fes = mesh_boundary(fes, []);
        bcl = fe_select(fens, bdry_fes, ...
            struct ('box',[a,a,-Inf,Inf,-Inf,Inf],'inflate',tolerance));
        traction.fes =subset(bdry_fes,bcl);;
        traction.traction= [0; 0; q0];
        traction.integration_rule =eltyd(eix).surface_integration_rule;
        model_data.boundary_conditions.traction{1} = traction;
        
        model_data.renumber=true;
        
        % Solve
        model_data =deformation_linear_statics(model_data);
        
        
        nc=[fenode_select(fens, struct('box', [a,a,0,b,0,0],...
            'inflate',tolerance))];
        bottom_edge_xyz=gather_values(model_data.geom,nc);
        bottom_edge_displacement=gather_values(model_data.u,nc);
        ns=[ns,model_data.u.nfreedofs];
        uzs=[uzs,bottom_edge_displacement(:,3)];
        uys=[uys,bottom_edge_displacement(:,2)];
        ys=[ys,bottom_edge_xyz(:,2)];
        
        if (graphics )
            figure
            clear options
            model_data.postprocessing.u_scale= u_scale;
            model_data.postprocessing.draw_mesh= false;
            %             model_data=deformation_plot_deformation(model_data);
            model_data.postprocessing.stress_component= 6;
            model_data.postprocessing.outputRm= eye(3);
            %model_data=deformation_plot_stress(model_data);
            model_data.postprocessing.stress_range= [-5800,5800];
            %             model_data=deformation_plot_stress(model_data);
            model_data=deformation_plot_stress_elementwise(model_data);
            anchor =[0,1.2*b,0]; tip = anchor +0.3*b*Rm(:, 1)';
            draw_arrow (model_data.postprocessing.gv,...
                anchor, tip-anchor, struct('facecolor', 'black'));
            anchor =[0,1.2*b,0]; tip = anchor +0.3*b*Rm(:, 2)';
            draw_arrow (model_data.postprocessing.gv,...
                anchor, tip-anchor, struct('facecolor', 'black'));
            draw_axes(model_data.postprocessing.gv,struct('length',(b/3), 'origin',[0,-b/2,0]));
            title(eltyd(eix).description)
            axis tight
            %             draw_polyline(model_data.postprocessing.gv,...
            %                 [anchor;  [tip(1),anchor(2),anchor(3)]; tip], [1,2;2,3],  struct('facecolor', 'black'))
            %             draw_polyline(model_data.postprocessing.gv,...
            %                 [anchor;  [anchor(1),tip(2),anchor(3)]; tip], [1,2;2,3],  struct('facecolor', 'black'))
            %             draw_polyline(model_data.postprocessing.gv,...
            %                 [anchor;  [anchor(1),anchor(2),tip(3)]; tip], [1,2;2,3],  struct('facecolor', 'black'))
        end
    end
    disp([']; %' eltyd(eix).description  ])
    disp(['hs=[' num2str(refs) '];' ])
    
    %         [xestim, beta, c, residual] = richextrapol(uzs(end-2:end),refs(end-2:end));
    if (~graphics )
        [~,ix]= sort(ys);
        plot(ys(ix),uzs(ix),name_to_style(eltyd(eix).description),'linewidth',3); hold on
        %                 plot(ns,uys,eltyd(eix).styl,'linewidth',3); hold on
        %                 plot(diff(ns),diff(uzs),eltyd(eix).styl,'linewidth',3); hold on
        figure (gcf); pause (1)
    end
    legends{end+1} =eltyd(eix).description;
end
if (~graphics )
legend (legends);
title (['Fiber at +/-45 degrees'])
labels(  'Distance along the bottom edge', 'Displacement in the direction of the load')
grid on
end
set_graphics_defaults
end
