% The initially twisted cantilever beam is one of the standard test
% problems for verifying the finite-element accuracy [1]. The beam is
% clamped at one end and loaded either with unit in-plane or 
% unit out-of-plane force at the other. The centroidal axis of the beam is
% straight at the undeformed  configuration, while its cross-sections are
% twisted about the centroidal axis from 0 at the clamped end to pi/2 at
% the free end. 
%
% Reference:
% Zupan D, Saje M (2004) On "A proposed standard set of problems to test
% finite element accuracy": the twisted beam. Finite Elements in Analysis
% and Design 40: 1445-1451.  
function  Twisted_beam
    % Parameters:
    E=0.29e8;
    nu=0.22;
    W=1.1;
    L=12;
    t= 0.32;
    nl=4; nt=1; nw=2;
    p=  1/W/t;
    %     Loading in the Z direction
    loadv=[0;0;p];dir=3;uzex=0.005424534868469; % Harder: 5.424e-3;
    %     Loading in the Y direction
    %         loadv=[0;p;0];dir=2;uzex=0.001753248285256; % Harder: 1.754e-3;
    graphics = ~true;
    u_scale=1000;
    
    
    prop=property_deformation_linear_iso(struct('E',E,'nu',nu));
    mater = material_deformation_linear_triax (struct('property',prop ));
    
    clear eltyd
    eix=1;
    
    
        eltyd(eix).description ='h8msgso';% tetrahedron
        eltyd(eix).mf =@H8_block;
        eltyd(eix).femmf =@(fes)femm_deformation_linear_h8msgso(struct('fes',fes,'material',mater,...
        'integration_rule',gauss_rule(struct('dim',3, 'order',2))));%;%tensprod_nq_rule(struct('dim',3, 'order',1))
        eltyd(eix).surface_integration_rule= gauss_rule(struct('dim',2, 'order',2));
        eltyd(eix).styl='b^-';
        eix=eix+1;
    
        %         eltyd(eix).description ='C8MS';% tetrahedron
        %         eltyd(eix).mf =@C8_block;
        %         eltyd(eix).femmf =@(fes)femm_deformation_nonlinear_c8ms(struct('fes',fes,'material',mater,...
        %         'integration_rule',tet_rule(struct('npts',1))));
        %         eltyd(eix).surface_integration_rule=tri_rule(struct('npts',1));
        %         eltyd(eix).styl='b^-';
        %         eix=eix+1;
    
    %     eltyd(eix).description ='H64';
    %     eltyd(eix).mf =@H64_block;
    %     eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
    %         'integration_rule',gauss_rule(struct('dim',3, 'order',4))));
    %     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',4));
    %     eltyd(eix).styl='k*--';
    %     eix=eix+1;
    
    %     eltyd(eix).description ='T10MS';% tetrahedron
    %     eltyd(eix).mf =@T10MS_block;
    %     eltyd(eix).femmf =@(fes)femm_deformation_nonlinear_t10ms(struct('fes',fes,'material',mater,'integration_rule',tet_rule(struct('npts',4))));
    %     eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
    %     eltyd(eix).styl='b^-';
    %     eix=eix+1;
    
%     eltyd(eix).description ='H27';
%     eltyd(eix).mf =@H27_block;
%     eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
%         'integration_rule',gauss_rule(struct('dim',3, 'order',3))));
%     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',4));
%     eltyd(eix).styl='ks--';
%     eix=eix+1;
    
    eltyd(eix).description ='T10';% tetrahedron
    eltyd(eix).mf =@T10_block;
    eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes,'material',mater,'integration_rule',tet_rule(struct('npts',4))));
    eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
    eltyd(eix).styl='k^--';
    eix=eix+1;
    
%     eltyd(eix).description ='T10-SRI';
%     eltyd(eix).mf =@T10_block;
%     eltyd(eix).femmf =@(fes)femm_deformation_linear_sri(struct('fes',fes, 'material',mater,...
%         'integration_rule_volumetric',tet_rule(struct('npts',1)),...
%         'integration_rule_deviatoric',tet_rule(struct('npts',4))));
%     eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
%     eltyd(eix).styl='kv-';
%     eix=eix+1;
    
    eltyd(eix).description ='H20R';
    eltyd(eix).mf =@H20_block;
    eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
        'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
    eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    eltyd(eix).styl='ro--';
    eix=eix+1;
    
    %         % Selective reduced integration hexahedron
    eltyd(eix).description ='H8-SRI';
    eltyd(eix).mf =@H8_block;
    eltyd(eix).femmf =@(fes)femm_deformation_linear_sri(struct('fes',fes, 'material',mater,...
        'integration_rule_volumetric',gauss_rule(struct('dim',3, 'order',1)),...
        'integration_rule_deviatoric',gauss_rule(struct('dim',3, 'order',2))));
    eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    eltyd(eix).styl='md--';
    eix=eix+1;
    
   mix = 1;
    mesd(mix).ref=1:4; %[1,2,4];
    mix = mix+1;
    
    for eix = 1:length(eltyd)
        ns=[];
        
        for mix = 1:length(mesd)
            uzs=[];
            nfreedofs = [];
            
            for     ref=mesd(mix).ref;
                %% Create the mesh and initialize the geometry
                [fens,fes]= eltyd(eix).mf(L,W,t, nl*ref,nw*ref,nt*ref);%max([round(nt*ref/2),2]));
                xy=fens.xyz;
                for i=1:count (fens)
                    a=xy(i,1)/L*(pi/2); y=xy(i,2)-(W/2); z=xy(i,3)-(t/2);
                    xy(i,:)=[xy(i,1),y*cos(a)-z*sin(a),y*sin(a)+z*cos(a)];
                end
                fens.xyz=xy;
                
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
                essential.node_list = fenode_select (fens,struct ('box',[0 0 -100*W 100*W -100*W 100*W],'inflate',0.001*t));
                model_data.boundary_conditions.essential{1} = essential;
                
                clear traction
                bdry_fes = mesh_boundary(fes, []);
                bcl = fe_select(fens, bdry_fes, ...
                    struct ('box',[L L -100*W 100*W -100*W 100*W],'inflate',0.0001*t));
                traction.fes =subset(bdry_fes,bcl);;
                traction.traction= loadv;
                traction.integration_rule =eltyd(eix).surface_integration_rule;
                model_data.boundary_conditions.traction{1} = traction;
                
                % Solve
                model_data =deformation_linear_statics(model_data);
                
                if (graphics )
                    clear options
                    options.u_scale= u_scale;
                    options=deformation_plot_deformation(model_data, options);
                end
                
                enl=fenode_select (fens,struct ('box',[L L -100*W 100*W -100*W 100*W],'inflate',0.01*t));
                uv=gather_values (model_data.u,enl);
                uz=max(uv(:,dir));
                disp (['%     Normalized displacement =' num2str(uz/uzex)])
                uzs =[uzs uz] ;
                nfreedofs = [nfreedofs,model_data.u.nfreedofs];
                ns=[ns,nl*ref] ;
            end
        end
        format long
        disp(['Data{end+1}=[']);
        disp([ns;uzs]);
        disp(['];' 'Style{end+1}=''' eltyd(eix).styl ''';' 'Legends{end+1} =''' eltyd(eix).description ''';']);
        plot(ns,uzs,eltyd(eix).styl);  hold on
    end
end