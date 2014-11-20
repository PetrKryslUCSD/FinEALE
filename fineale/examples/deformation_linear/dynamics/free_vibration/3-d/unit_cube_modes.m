%
% Vibration modes of unit cube  of almost incompressible material.
%
% Reference: Puso MA, Solberg J (2006) A stabilized nodally integrated
% tetrahedral. International Journal for Numerical Methods in
% Engineering 67: 841-867.
function unit_cube_modes
    pu=physical_units_struct;
    % Parameters:
    E = 1*pu.PA;% 210e3 MPa
    nu = 0.499;
    rho= 1*pu.KG/pu.M^3;
    a=1*pu.M; b=a; h= a;
    % Mesh refinement
    n1=1;% How many element edges per side?
    na= n1; nb= n1; nh =n1;
    %     Shall we produce graphics?
    graphics = ~false;
    make_frames = false;% Shall we save the frames as images?
    
    
    prop=property_deformation_linear_iso(struct('E',E,'nu',nu,'rho',rho));
    mater = material_deformation_linear_triax (struct('property',prop ));
    
    clear eltyd
    eix=1;
    
    
    %             eltyd(eix).description ='T10';% tetrahedron
    %             eltyd(eix).mf =@T10_block;
    %             eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes,'material',mater,...
    %                 'integration_rule',tet_rule(struct('npts',4))));
    %             eltyd(eix).styl='gd-';
    %             eix=eix+1;
    %
    eltyd(eix).description ='T10-SRI';
    eltyd(eix).mf =@T10_block;
    eltyd(eix).femmf =@(fes)femm_deformation_linear_sri(struct('fes',fes, 'material',mater,...
        'integration_rule_volumetric',tet_rule(struct('npts',1)),...
        'integration_rule_deviatoric',tet_rule(struct('npts',4)),...
        'integration_rule',tet_rule(struct('npts',4))));
    eltyd(eix).styl='gd-';
    eix=eix+1;
    
    %     % Selective reduced integration hexahedron
    %     eltyd(eix).description ='H8-SRI';
    %     eltyd(eix).mf =@H8_block;
    %     eltyd(eix).femmf =@(fes)model_deformation_linear_sri(struct('fes',fes, 'material',mater,...
    %         'integration_rule_volumetric',gauss_rule(struct('dim',3, 'order',1)),...
    %         'integration_rule_deviatoric',gauss_rule(struct('dim',3, 'order',2))));
    %      eltyd(eix).styl='md-';
    %     eix=eix+1;
    
    %     eltyd(eix).description ='H20R';
    %     eltyd(eix).mf =@H20_block;
    %     eltyd(eix).femmf =@(fes)femm_deformation_linear(struct('fes',fes, 'material',mater,...
    %         'integration_rule',gauss_rule(struct('dim',3, 'order',2))));
    %     eltyd(eix).surface_integration_rule=gauss_rule(struct('dim',2, 'order',2));
    %     eltyd(eix).styl='ro-';
    %     eix=eix+1;
    
    
    
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
        model_data.omega_shift=(2*pi*0.01) ^ 2;
        model_data.use_factorization= true;
        % Solve
        model_data = deformation_linear_modal_analysis(model_data);
        
        f=model_data.Omega'/2/pi;
        disp ([eltyd(eix).description  ': '])
        disp(['Frequency ' num2str(abs(f)) ' [Hz]' ]);
        disp(['Frequency ' num2str((f)) ' [Hz]' ]);
        %         disp(['             f/f_analytical=' num2str(f./f_fundamental(7:10)*100) '%']);
        %         clf;
        
        model_data.postprocessing.u_scale= 2;
        model_data.postprocessing.modelist= 1:13;
        model_data.postprocessing.save_frame= make_frames;
        model_data.postprocessing.frame_name= [eltyd(eix).description '-n1-' num2str(n1)];
        model_data=deformation_plot_modes(model_data);
        
        
    end
end