% Vibration analysis of laminated composite plate
% Reference:
% Analytical solutions for free vibration of laminated composite and sandwich
% plates based on the higher-order refined theory, T. Kant and K. Swaminathan,
% Composite Structures 53  (2001) 73-85
function KS_sandwich_vibration_ss
    
    % Material: The data from the above cited article. The ratio of the Young's
    % moduli is variable, and the number of layers also varies.
    % E1/E2 = 3, 10, 40
    
%     % One particular case from the paper:
%     L_t =5;
%     E1_E2 = 3;% this ratio is variable
%     % Anti-symmetric layup
%     angles =[0 90 0 90];
%     omega_bar=6.5455; % the fundamental frequency multiplier
    
%         % One particular case from the paper:
%         L_t =10;
%         E1_E2 = 40;% this ratio is variable
%         % Anti-symmetric layup
%         angles =[0 90 ];
%         omega_bar=10.4319; % the fundamental frequency multiplier
    
        %         % One particular case from the paper:
        %         L_t =20;
        %         E1_E2 = 40;% this ratio is variable
        %         % Anti-symmetric layup
        %         angles =[0 90 ];
        %         omega_bar=11.0663; % the fundamental frequency multiplier
    
        % One particular case from the paper:
        L_t =100;
        E1_E2 = 40;% this ratio is variable
        % Anti-symmetric layup
        angles =[0 90 ];
        omega_bar=11.2988; % the fundamental frequency multiplier
 
    % More anisotropy:
    % L_t =5;
    % E1_E2 = 40;% this ratio is variable
    % % Anti-symmetric layup
    % angles =[0 90 0 90];
    % omega_bar=10.6798; % the fundamental frequency multiplier
    
    % Increased number of layers:
    %     L_t =5;
    %     E1_E2 = 40;% this ratio is variable
    %     % Anti-symmetric layup
    %     angles =[0 90 0 90  0 90 0 90 0 90];
    %     omega_bar=11.6245; % the fundamental frequency multiplier
    
    % Compute  all the other parameters
    E1=181000; E2=E1/E1_E2; E3=E2; G12=0.6*E2;  G13=0.6*E2; G23=0.5*E2;% MPa
    nu12= 0.25; nu13= 0.25; nu23= 0.25; rho =1.5e-9;
    L=200; W=200; t=L/L_t;% mm
    
    % Fundamental frequency
    omega_fundamental =omega_bar*t/W^2*sqrt(E2/rho);
    f_fundamental=omega_fundamental/2/pi;
    disp(['Analytical f_fundamental =',num2str(f_fundamental)])
    
    nLayers =length(angles);
    [nL,nW] =adeal(2*[4,4]);
    nts= 1*ones(length(angles));% number of elements per layer
    ts= t/length(angles)*ones(length(angles));% layer thicknesses
    graphics = ~false;
    
    
    
    function Rm = LayerRm(XYZ, ts, label)% label equals the layer number here
        Rm= rotmat(angles(label)/180*pi* [0,0,1]);
    end
    % Mesh
    [fens,fes] = H8_composite_plate(L,W,ts,nL,nW,nts);;
    [fens,fes] = H8_to_H20(fens,fes);
    %     gv=drawmesh( {fens,subset(fes,fe_select(fens,fes,struct('label', 1)))},'fes', 'facecolor','r');
    %     gv=drawmesh( {fens,subset(fes,fe_select(fens,fes,struct('label', 2)))},'gv',gv,'fes', 'facecolor','g');
    %     gv=drawmesh( {fens,subset(fes,fe_select(fens,fes,struct('label', 3)))},'gv',gv,'fes', 'facecolor','b');
    %     gv=drawmesh( {fens,subset(fes,fe_select(fens,fes,struct('label', 4)))},'gv',gv,'fes', 'facecolor','m');
    
    
    % Compose the model data
    clear model_data
    model_data.fens =fens;
    
    clear region
    region.property = 'orthotropic';
    region.rho =rho;
    region.E1 =E1;
    region.E2 =E2;
    region.E3 =E3;
    region.G12=G12;
    region.G13=G13;
    region.G23=G23;
    region.nu12=nu12;
    region.nu13=nu13;
    region.nu23=nu23;
    region.fes= fes;
    region.integration_rule = gauss_rule (struct('dim', 3, 'order', 2));
    region.Rm =@LayerRm;
    model_data.region{1} =region;
    
    clear essential
    essential.component= [2,3];
    essential.fixed_value= 0;
    essential.node_list = [fenode_select(fens, struct('box', [0,0,-Inf,Inf,-Inf,Inf],...
        'inflate',0.001*t)),fenode_select(fens, struct('box', [L,L,-Inf,Inf,-Inf,Inf],...
        'inflate',0.001*t))];
    model_data.boundary_conditions.essential{1} = essential;
    
    clear essential
    essential.component= [1,3];
    essential.fixed_value= 0;
    essential.node_list = [fenode_select(fens, struct('box', [0,L,0,0,-Inf,Inf],...
        'inflate',0.001*t)),fenode_select(fens, struct('box', [0,L,W,W,-Inf,Inf],...
        'inflate',0.001*t))];
    model_data.boundary_conditions.essential{2} = essential;
    
    model_data.neigvs=6;
    model_data.omega_shift=(2*pi*40) ^ 2;
    model_data.use_factorization= true;
    % Solve
    tic; model_data = deformation_linear_modal_analysis(model_data);
    toc
    
    f=model_data.Omega(1)/2/pi;
    disp(['  Eigenvector ' num2str(1) ' frequency ' num2str(f) ' [Hz]' ', f/f_analytical=' num2str(f/f_fundamental)]);
    %         clf;
    
    clear options
    model_data.postprocessing.u_scale= 2;
    model_data.postprocessing.modelist= 1:model_data.neigvs;
    model_data=deformation_plot_modes(model_data);
    
    
end