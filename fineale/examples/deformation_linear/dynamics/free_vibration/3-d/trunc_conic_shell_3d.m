% Cylindrical  free-floating shell.  Free vibration modes.
%    See the analytical solution in the enclosed function.
%
%
function trunc_conic_shell_3d
    u=  physical_units_struct;
    E=205000*u.MEGA*u.PA;% Young's modulus
    nu=0.3;% Poisson ratio
    rho=7850*u.KG*u.M^-3;% mass density
    OmegaShift=(2*pi*100) ^ 2;
    h=0.05*u.M;
    l=10*h;
    Rm=h/0.2;
    % Cylinder
    psi  =0;
    nh=2; nl =4; nc=8;
    tolerance=h/nh/100;
    neigvs=20;
    
    graphics = ~false;
    u_scale=1;
    
    % Mesh
    [fens,fes]=H8_block(h,l,2*pi,nh,nl,nc);
    [fens,fes] = H8_to_H64(fens,fes);
    xyz=fens.xyz;
    for i=1:count (fens)
        x=xyz(i,1); y=xyz(i,2); z=xyz(i,3);
        xyz(i,:)= [x+Rm-h/2, y-l/2, 0];
        xyz(i,:) = xyz(i,:)*[cos(psi*pi/180),sin(psi*pi/180),0;-sin(psi*pi/180),cos(psi*pi/180),0;0,0,1]*rotmat([0, z, 0]);
    end
    fens.xyz= xyz;
    [fens,fes] = merge_nodes(fens, fes,  tolerance);
    
    
    % Compose the model data
    clear model_data
    model_data.neigvs =20;
    model_data.omega_shift =200*2*pi;
    
    model_data.fens =fens;
    
    clear region
    region.E =E;
    region.nu =nu;
    region.rho =rho;
    region.fes= fes;
    region.integration_rule = gauss_rule (struct('dim', 3, 'order', 4));
    region.Rm =[];
    model_data.region{1} =region;
    
    
    % Solve
    analytical_solution
    model_data = deformation_linear_modal_analysis(model_data);
    
    model_data.postprocessing.u_scale= u_scale;
    model_data.postprocessing.modelist= 7;
    model_data=deformation_plot_modes(model_data);
    
    function  analytical_solution
        
        format short
        disp( ['Analytical natural frequencies in Hz from LADEFOGED (1988)'])
        c= [0.1496, 0.1942, 0.7564, 1.0553];
        omega=c/Rm*sqrt(E/(1-nu^2)/rho);
        f=omega/2/pi;
        num2str(f,5)
        
        
        disp( ['Analytical natural frequencies in Hz from GLADWELL AND VIJAY (1975)'])
        H= 0.2; L= 2.0;
        c= [0.2529 %first symmetric frequency for circumferential wave number = 2
            0.3282 %first anti-symmetric frequency for circumferential wave number = 2
            1.2592 % first anti-symmetric frequency for circumferential wave number = 1
            1.2786 % second symmetric frequency for circumferential wave number = 2
            1.5765 % first symmetric frequency for circumferential wave number = 0
            1.6052 % first anti-symmetric frequency for circumferential wave number = 0
            ];
        omega=c'/Rm*sqrt(E/2/(1+nu)/rho);
        f=omega/2/pi;
        num2str(f,5)
    end
end
