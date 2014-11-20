% Taylor clamped beam with end shear force
%
function rltb_3D
    E=1000;
    nu=0.4999;uzex =-12.0935378981478;%12.2334;%Taylor data % from H64
    %         nu=0.29;uzex=-12.2506307034569;% 12.3939;%Testing
        %             nu=0.49999999;uzex=- 12.2450;%%Testing
%     nu=0.4999999; uzex=-12.0805262692106; 
    %         nu=0.49;%Testing
    W=2.5;
    H=5;
    L= 50;
    htol=min([L,H,W])/1000;
    magn = -0.2*12.2334/4;
    Force =magn*W*H*2;
    Force*L^3/(3*E*W*H^3*2/12);
    graphics = ~true;
    stress_range = [-30,30];
    u_scale=1;
    
    clc;
    disp(['% ' mfilename]);
    disp(['nu=' num2str(nu,12) ';']);
    disp(['uzex=' num2str(uzex,12) ';']);
    
    
    prop=property_deformation_linear_iso(struct('E',E,'nu',nu));
    mater = material_deformation_linear_triax (struct('property',prop ));
                
    mix = 1;
    mesd(mix).mult=2;
    mesd(mix).ns=[2,4,8];
    mesd(mix).ns=[1,2,4];
    mix = mix+1;
    
    clear eltyd
    eix=1;
 
    eltyd(eix).description ='T10';% tetrahedron
    eltyd(eix).mf =@T10_block;
    eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes,'material',mater,'integration_rule',tet_rule(struct('npts',4))));
    eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
    eltyd(eix).styl='k^--';
    eix=eix+1;
    
    eltyd(eix).description ='H8-SRI';% tetrahedron
    eltyd(eix).mf =@H8_block;
    eltyd(eix).femmf =@(fes)femm_deformation_linear_sri(struct('fes',fes,'material',mater,...
        'integration_rule_volumetric',gauss_rule(struct('dim',3, 'order',1)),...
        'integration_rule_deviatoric',gauss_rule(struct('dim',3, 'order',2))));
    eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    eltyd(eix).styl='gs--';
    eix=eix+1;
    
    
    eltyd(eix).description ='H20R';
    eltyd(eix).mf =@H20_block;
    eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
        'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
    eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    eltyd(eix).styl='ro--';
    eix=eix+1;
    
    eltyd(eix).description ='H20';
    eltyd(eix).mf =@H20_block;
    eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
        'integration_rule',gauss_rule(struct('dim',3, 'order',3))));
    eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    eltyd(eix).styl='r.--';
    eix=eix+1;
    
    legends ={};
    for eix = 1:length(eltyd)
        
        for mix = 1
            ns=[];      nhs=[];      uzs=[];
            for     n=mesd(mix).ns;
                mult =mesd(mix).mult;
                nW=n; nL =mult*n; nH=2*n;
            %% Create the mesh and initialize the geometry
                [fens,fes]= eltyd(eix).mf(W,L,H,nW,nL,nH);
                %                     drawmesh({fens,fes},'shrink', 0.8,'facecolor','red');
                nhs=[nhs,nH];
                
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
                essential.node_list = fenode_select(fens,struct('box',[0 W 0 0 0 H],'inflate',htol));
                model_data.boundary_conditions.essential{1} = essential;
                
                clear essential
                essential.component= [1];
                essential.fixed_value= 0;
                essential.node_list = fenode_select(fens,struct('box',[W W 0 L 0 H],'inflate',htol));
                model_data.boundary_conditions.essential{2} = essential;
                
                clear traction
                bdry_fes = mesh_boundary(fes, []);
                bcl = fe_select(fens, bdry_fes, ...
                    struct ('box',[0 W L L 0 H],'inflate',htol));
                traction.fes =subset(bdry_fes,bcl);;
                traction.traction= [0;0;magn];
                traction.integration_rule =eltyd(eix).surface_integration_rule;
                model_data.boundary_conditions.traction{1} = traction;
                
                % Solve
                model_data =deformation_linear_statics(model_data);
               
                
                ebc_fenids=fenode_select(fens,struct('box',[0 W L L 0 H]));
                uv=gather_values(model_data.u,ebc_fenids);
                uz=sum(uv(:, 3))/length(ebc_fenids);
                disp (['% ' eltyd(eix).description  ': uz =' num2str(uz) ', ' num2str(uz/uzex*100) ' %'])
                ns=[ns,model_data.u.nfreedofs];
                uzs =[uzs uz];
                
                if (graphics )
                    clear options
                    model_data.postprocessing.u_scale= u_scale;
                    model_data.postprocessing.stress_component=2;
%                     options=deformation_plot_stress(model_data, options)
                    model_data=deformation_plot_deformation(model_data);
                else
                    %                     plot(ns,uzs/uzex,eltyd(eix).styl,'linewidth',3); hold on
                    %                     figure (gcf); pause (1)
                end
                
            end
        end
        if (graphics )
        else
            %             semilogx(ns,uzs/uzex,eltyd(eix).styl,'linewidth',3); hold on
            %             loglog(ns,abs((uzs-uzex)/uzex),eltyd(eix).styl,'linewidth',3); hold on
            %             title(['Taylor, \nu=' num2str(nu,12), ', mult =' [ num2str( [mult])']'])
             plot(nhs,uzs/uzex,eltyd(eix).styl,'linewidth',3); hold on
            figure (gcf); pause (1)
        end
        legends{end+1} =eltyd(eix).description;
        format long
        disp(['Data{end+1}=[']);
        disp([num2str(ns,15)]);
        disp([num2str(uzs,15)]);
        disp(['];' 'Style{end+1}=''' eltyd(eix).styl ''';' 'Legends{end+1} =''' eltyd(eix).description ''';']);
        [xestim, beta, c, residual] = richextrapol(abs(uzs),mesd(mix).ns) ;
        disp(['% nu=' num2str(nu,15) '; ' 'uzex=-' num2str(xestim,15) '; ']);
    end
    legend (legends);
    labels('Number of elements per height',  'Normalized deflection')
    grid on
    set_graphics_defaults
end