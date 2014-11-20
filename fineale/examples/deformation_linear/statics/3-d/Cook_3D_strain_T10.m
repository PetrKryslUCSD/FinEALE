% Cook's membrane under plane strain conditions, 3-D version, tetrahedral elements
function  Cook_3D_strain_T10
    % Parameters:
    lambda = 0.75e4;
    %     lambda = 0.75e9;
    G= 0.375;
    convutip=16.432584376299037;% NICE-T6 t3block2d with Richardson
    nu =1/2*lambda/(G+lambda);
    E=G*(3*lambda+2*G)/(G+lambda);
    tol=48/10000;
    thickness=48/100;
    magn=1;
    aspect=1;
    graphics = ~true;
    u_scale=1;
    
    ns=[1,2,4,8,16,32];
    ns=[1,2,4,8,16];
    
    format long
    clc
    disp(['% ' mfilename]);
    disp(['nu=' num2str(nu,12)]);
    
    prop=property_deformation_linear_iso(struct('E',E,'nu',nu));
    mater = material_deformation_linear_triax (struct('property',prop ));
    
    clear eltyd
    eix=1;
    
    
    eltyd(eix).description ='T10';% tetrahedron
    eltyd(eix).mf =@T10_block;
    eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes,'material',mater,'integration_rule',tet_rule(struct('npts',4))));
    eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
    eltyd(eix).styl='k^--';
    eix=eix+1;
    
    eltyd(eix).description ='T10-SRI';
    eltyd(eix).mf =@T10_block;
    eltyd(eix).femmf =@(fes)femm_deformation_linear_sri(struct('fes',fes, 'material',mater,...
        'integration_rule_volumetric',tet_rule(struct('npts',1)),...
        'integration_rule_deviatoric',tet_rule(struct('npts',4))));
    eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
    eltyd(eix).styl='kv-';
    eix=eix+1;
    
    
    for eix = 1:length(eltyd)
        uzs=[];
        nfreedofs = [];
        
        for n=ns  %,16,24,32,48
            % Create the mesh and initialize the geometry
            [fens,fes]= eltyd(eix).mf(48,44,thickness/n, n, aspect*n,1);%max([round(nt*ref/2),2]));
            dxy=min(48,44)/n/aspect/100;
            xy=fens.xyz;
            for i=1:count (fens)
                xy(i,:)=[xy(i,1),xy(i,2) + (xy(i,1)/48)*(44 -xy(i,2)/44*28),xy(i,3)];
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
            essential.node_list = fenode_select (fens,struct ('box',[0,0,0, 44,0,thickness/n],'inflate',tol));
            model_data.boundary_conditions.essential{1} = essential;
            
            clear essential
            essential.component= [3];
            essential.fixed_value= 0;
            essential.node_list = [fenode_select(fens,struct('box',[-inf,inf,-inf,inf,0,0],'inflate',tol)),...
                fenode_select(fens,struct('box',[-inf,inf,-inf,inf,thickness/n,thickness/n],'inflate',tol))];
            model_data.boundary_conditions.essential{2} = essential;
            
            clear traction
            bdry_fes = mesh_boundary(fes, []);
            bcl = fe_select(fens, bdry_fes, ...
                struct ('box',[48,48,44,60,0,thickness/n],'inflate',tol));
            traction.fes =subset(bdry_fes,bcl);;
            traction.traction= [0;magn/(60-44);0];
            traction.integration_rule =eltyd(eix).surface_integration_rule;
            model_data.boundary_conditions.traction{1} = traction;
            
            % Solve
            model_data =deformation_linear_statics(model_data);
            
            if (graphics )
                clear options
                options.u_scale= u_scale;
                options=deformation_plot_deformation(model_data, options);
            end
            
            enl=fenode_select (fens,struct ('box',[48,48,52,52,0,thickness],'inflate',tol));
            uv=gather_values (model_data.u,enl);
            utip=mean(uv(:,2));
            disp (['%     Normalized displacement =' num2str(utip/convutip)])
            uzs =[uzs utip] ;
            nfreedofs = [nfreedofs,model_data.u.nfreedofs];
        end
        format long
        disp(['Data{end+1}=[']);
        disp([ns;uzs]);
        disp(['];' 'Style{end+1}=''' eltyd(eix).styl ''';' 'Legends{end+1} =''' eltyd(eix).description ''';']);
        plot(ns,uzs,eltyd(eix).styl);  hold on
    end
end