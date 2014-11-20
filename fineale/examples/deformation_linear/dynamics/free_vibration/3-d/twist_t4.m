% Transverse isotropic, twisted fiber rod. T4 elements.
function twist_t4
    nR = [2, 3, 4];
    nt = [30, 40, 50, 60];
    fs =zeros(length(nR),length(nt));
    rand('state',0);% try to comment out this line and compare
    %                   results for several subsequent runs
    for i=1:length(nR)
        for j=1:length(nt)
            fs(i,j) =do_twist_t4(nR(i),nt(j))
        end
    end
end


function frequency1 =do_twist_t4(nR,nt)
    % graphite polymer composite
    E1=155000;
    E2=12100;
    G12=4400;
    nu12=0.248;
    nu23=0.458;
    rho=5e-9;
    integration_order=1;
    R= 2.5;
    t= 100.0;
    twist_angle= 15/360*2*pi;
    graphics = false;
    
    function Rm = default (XYZ, ts, label)
        Rm= eye(3);
    end
    function Rm = twist (XYZ, ts,  label)
        r= norm(XYZ( 2:3));
        if r>0
            y=XYZ(2);z=XYZ(3);
            e2 = [0, y/r, z/r]';
            e3 = skewmat([1, 0, 0])*e2;
            Rm= [[1, 0, 0]',e2,e3];
            Rm=rotmat(r/R*twist_angle*skewmat(e2))*Rm;
        else
            Rm= eye(3);
        end
    end
    % Mesh
    [fens,fes] = T4_cylinderdel(t,R, nt,nR);
    % Material
    prop = property_deformation_linear_transv_iso (...
        struct('E1',E1,'E2',E2,'G12',G12,'nu12',nu12,'nu23',nu23,'rho',rho));
    mater = material_deformation_linear_triax (struct('property',prop ));
    %     Make the finite element model machine
    femm = femm_deformation_linear (struct ('material',mater, 'fes',fes,...
        'integration_rule',tet_rule (struct('npts',1)),'Rm',@twist));
    
    % Geometry
    geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
    % Define the displacement field
    u   = 0*geom; % zero out
    % Apply EBC's
    nl=[fenode_select(fens, struct('box', [0,0,-R,R,-R,R],...
        'inflate',0.001)),...
        fenode_select(fens, struct('box', [t,t,-R,R,-R,R],...
        'inflate',0.001))];
    u   = set_ebc(u,nl,nl*0+1, [], nl*0);
    u   = apply_ebc (u);
    % Number equations
    u   = numberdofs (u);
    % Assemble the system matrix
    K = stiffness(femm, sysmat_assembler_sparse, geom, u);
    M = mass(femm, sysmat_assembler_sparse, geom, u);
    u.nfreedofs
    
    %
    neigvs = 20;
    options.issym= true;
    options.tol=1e-80;
    options.maxit= 500;
    options.disp= 0;
    [W,Omega]=eigs(K,M,neigvs,'SM', options);
    frequency1=sqrt(Omega(1, 1))/2/pi;
    
    % Plot
   if graphics
        bg=mesh_boundary(fes);
        gv=graphic_viewer;
        gv=reset (gv,[]);
        scale=125;
        w=clone(u,'w'); % make a copy of u
        for i=1:1
            disp(['  Eigenvector ' num2str(i) ' frequency ' num2str(sqrt(Omega(i,i))/2/pi) ]);
            clf;
            w = scatter_sysvec(w, W(:,i));
            wv = magnitude(w);
            wmag = wv.values;
            dcm=data_colormap(struct ('range',[min(wmag),max(wmag)], 'colormap',jet));
            colors=map_data(dcm, wmag);
            colorfield = nodal_field(struct ('name',['cf'], 'dim', 3, 'data',colors));
            gv=reset (gv,[]);
            set(gca,'FontSize',16)
            camset(gv,[885.2912  -77.5338  228.2605   46.4468    3.2376    2.4532   -0.2576    0.0248    0.9659    2.6098]);
            draw(bg,gv, struct ('x', geom,'u',+scale*w,'colorfield',colorfield));
            draw_axes(gv)
            %         text(3.1*R,4.1*R,4.1*R,['\omega_' num2str(i) '=' num2str(sqrt (Omega(i,i))) ],'FontSize',24);
            axis off
            %     saveas(gcf, ['twist_t4-' num2str(i) '.png'], 'png');
            pause(2); % next eigenvector
        end
    end
end
