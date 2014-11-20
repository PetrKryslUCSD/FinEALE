% Cook's membrane under plane strain conditions, 3-D version
function  Cook_3D_strain
% Parameters:
lambda = 0.75e4;
G= 0.375;
convutip=16.432584376299037;% NICE-T6 t3block2d with Richardson
nu =1/2*lambda/(G+lambda);
E=G*(3*lambda+2*G)/(G+lambda);
tol=48/10000;
thickness=48/1;
magn=1;
aspect=1;
graphics = ~true;
Plot_stress = true;
u_scale=1;

ns=[1,2,4,8,16];

format long
clc
disp(['% ' mfilename]);
disp(['nu=' num2str(nu,12) '; ']);

prop=property_deformation_linear_iso(struct('E',E,'nu',nu));
mater = material_deformation_linear_triax (struct('property',prop ));

clear eltyd
eix=1;


eltyd(eix).description ='H8-Bbar';% tetrahedron
eltyd(eix).mf =@H8_block;
eltyd(eix).femmf =@(fes)femm_deformation_linear_bbar(struct('fes',fes,'material',mater,...
    'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',1)),...
    'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',2))));
eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
eix=eix+1;

% eltyd(eix).description ='H8-GSRI';% tetrahedron
% eltyd(eix).mf =@H8_block;
% eltyd(eix).femmf =@(fes)femm_deformation_linear_gsri(struct('fes',fes,'material',mater,...
%     'integration_rule_constrained',gauss_rule(struct('dim',3, 'order',1)),...
%     'integration_rule_unconstrained',gauss_rule(struct('dim',3, 'order',2))));
% eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
% eix=eix+1;

% eltyd(eix).description ='T10';% tetrahedron
% eltyd(eix).mf =@T10_block;
% eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes,'material',mater,'integration_rule',tet_rule(struct('npts',4))));
% eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
% eltyd(eix).styl='k^--';
% eix=eix+1;

%         eltyd(eix).description ='T10-SRI';
%         eltyd(eix).mf =@T10_block;
%         eltyd(eix).femmf =@(fes)femm_deformation_linear_sri(struct('fes',fes, 'material',mater,...
%             'integration_rule_volumetric',tet_rule(struct('npts',1)),...
%             'integration_rule_deviatoric',tet_rule(struct('npts',4))));
%         eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
%         eltyd(eix).styl='kv-';
%         eix=eix+1;


%             eltyd(eix).description ='H20R';
%             eltyd(eix).mf =@H20_block;
%             eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
%                 'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
%             eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
%             eltyd(eix).styl='ro--';
%             eix=eix+1;
%
%         eltyd(eix).description ='H20';
%         eltyd(eix).mf =@H20_block;
%         eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
%             'integration_rule',gauss_rule(struct('dim',3, 'order',3))));
%         eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
%         eltyd(eix).styl='r.--';
%         eix=eix+1;
%
%     eltyd(eix).description ='H27';
%     eltyd(eix).mf =@H27_block;
%     eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
%         'integration_rule',gauss_rule(struct('dim',3, 'order',3))));
%     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
%     eltyd(eix).styl='rh--';
%     eix=eix+1;
%


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
            if ( Plot_stress )
                model_data.postprocessing.u_scale= u_scale;
                model_data.postprocessing.stress_component= 'pressure';
                model_data.postprocessing.use_spr = ~true;;
                model_data.postprocessing.stress_range  =[-0.2, 0.2];
                model_data=deformation_plot_stress(model_data);
            else
                model_data.postprocessing.u_scale= u_scale;
                model_data=deformation_plot_deformation(model_data);
            end
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
    disp(['];' 'Style{end+1}=''' name_to_style(eltyd(eix).description) ''';' 'Legends{end+1} =''' eltyd(eix).description ''';']);
    if ( ~graphics )
        plot(ns,uzs,name_to_style(eltyd(eix).description));  hold on
    end
end
end