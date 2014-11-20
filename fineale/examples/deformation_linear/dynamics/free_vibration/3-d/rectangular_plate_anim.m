% Vibration of a rectangular plate.
% Reference: C.M. Wang, W.X. Wu, C. Shu, T. Utsunomiya, LSFD method
% for accurate vibration modes and modal stress-resultants of freely
% vibrating plates that model VLFS, Computers and Structures 84 (2006)
% 2329–2339
function rectangular_plate
    
    pu=physical_units_struct;
    % Parameters:
    E = 210e3*pu.MEGA*pu.PA;% 210e3 MPa
    nu = 0.3;
    rho= 7850*pu.KG/pu.M^3;
    a=4.00*pu.M; b=1.00*pu.M; h= 0.1*pu.M;
    % Mesh refinement
    na= 13; nb=  9; nh =3;
    na= 29; nb=  9; nh =3;
    na= 10; nb=  5; nh =2;
    %     na= 13; nb=  4; nh =1;
    % na= 5; nb=  4; nh =2;
    
    % Fundamental frequency
    f_fundamental =[0 0 0  0 0 0 33.78, 82.28,  92.99, 170.06];
    graphics = ~false;
    
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
%         eltyd(eix).description ='T10';% tetrahedron
%         eltyd(eix).mf =@T10_block;
%         eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes,'material',mater,'integration_rule',tet_rule(struct('npts',4))));
%         eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
%         eltyd(eix).styl='k^--';
%         eix=eix+1;
    %
        eltyd(eix).description ='T10-SRI';
        eltyd(eix).mf =@T10_block;
        eltyd(eix).femmf =@(fes)femm_deformation_linear_sri(struct('fes',fes, 'material',mater,...
            'integration_rule_volumetric',tet_rule(struct('npts',1)),...
            'integration_rule_deviatoric',tet_rule(struct('npts',4)),...
            'integration_rule',tet_rule(struct('npts',4))));
        eltyd(eix).surface_integration_rule=tri_rule(struct('npts',3));
        eltyd(eix).styl='kv-';
        eix=eix+1;
        
    %
    %     eltyd(eix).description ='H20R';
    %     eltyd(eix).mf =@H20_block;
    %     eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
    %         'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
    %     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    %     eltyd(eix).styl='ro--';
    %     eix=eix+1;
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
    
    
    for eix = 1:length(eltyd)
        
        ns=[];
        uzs=[];
        
        %% Create the mesh and initialize the geometry
        [fens,fes]= eltyd(eix).mf(a,b,h,na,nb,nh);
        %                     drawmesh({fens,fes},'shrink', 0.8,'facecolor','red');
        
        % Compose the model data
        clear model_data
        model_data.fens =fens;
        
        clear region
        region.fes= fes;
        region.femm= eltyd(eix).femmf(fes);
        model_data.region{1} =region;
        
        model_data.neigvs= 20;
        model_data.omega_shift=(2*pi*20) ^ 2;
        model_data.use_factorization= true;
        % Solve
        model_data = deformation_linear_modal_analysis(model_data);
        
        f=model_data.Omega(7:10)'/2/pi;
        disp ([eltyd(eix).description  ': '])
        disp(['  Eigenvector ' num2str(1) ' frequency ' num2str(f) ' [Hz]' ]);
        disp(['             f/f_analytical=' num2str(f./f_fundamental(7:10)*100) '%']);
      
        model_data.postprocessing.u_scale= 2;
        model_data.postprocessing.modelist= 7;
        model_data.postprocessing.animate= true;
        model_data.postprocessing.camera = [-16.3474  -23.6183   17.5886    2.0624    0.3739    0.1287         0         0    1.0000    4.0316];
        model_data=deformation_plot_modes(model_data);
        
        
    end
end