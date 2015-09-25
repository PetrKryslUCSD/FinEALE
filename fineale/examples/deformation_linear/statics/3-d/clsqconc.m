function clsqconc
  
% Clamped square plate with central force.
% !  Timoshenko: clamped square plate with center load
% ! w_max = 0.0056 * P * a^2 / (E * h^3 / (12 * (1 - nu^2))) =
% !     2.44608 for P=10, a=2, E=1e9, h=0.001, nu=0.3
% 
  disp('% Timoshenko: clamped square plate with center load');
% Parameters:
graphics =~false; % graphic output
% Mesh
neqs=[]; normalized_deflections = [];
for nt=2
    for na =4:4:20
        P=2; a=2; E=1e9; h=0.05; nu=0.3;
        w_max = 0.0056 * P * a^2 / (E * h^3 / (12 * (1 - nu^2)))
        scale = 0.5/w_max;
        % tic;
        [fens,fes] = H8_block (a/2,a/2,h,na,na,nt);
                [fens,fes] = H8_to_H20(fens,fes);
        ir=gauss_rule(struct('dim',3,'order',2));
        % [fens,fes] = t4block (a/2,a/2,h,na,na,nt);
        %         [fens,fes] = T4_to_T10(fens,fes);
        % ir= tet_rule(4);
        % Material
        prop = property_deformation_linear_iso (struct('E',E,'nu',nu));
        mater = material_deformation_linear_triax (struct('property',prop));
        % Finite element block
        femm = femm_deformation_linear(struct('material',mater,'fes',fes,...
            'integration_rule',ir));
        % Geometry
        geom = nodal_field(struct ('name',['geom'], 'dim', 3, 'fens',fens));
        % Define the displacement field
        u   = 0*geom; % zero out
        % Apply EBC's
        ebc_fenids=fenode_select (fens,struct('box',[0 0 0 a/2 0 h],'inflate',h/100));
        ebc_fixed=ebc_fenids*0+1;
        ebc_comp=ebc_fenids*0+1;
        ebc_val=ebc_fenids*0;
        u   = set_ebc(u, ebc_fenids, ebc_fixed, ebc_comp, ebc_val);
        ebc_fenids=fenode_select (fens,struct('box',[0 a/2 0 0 0 h],'inflate',h/100));
        ebc_fixed=ebc_fenids*0+1;
        ebc_comp=ebc_fenids*0+2;
        ebc_val=ebc_fenids*0;
        u   = set_ebc(u, ebc_fenids, ebc_fixed, ebc_comp, ebc_val);
        ebc_fenids=fenode_select (fens,struct('box',[0 a/2 a/2 a/2 0 h],'inflate',h/100));
        ebc_fixed=ebc_fenids*0+1;
        ebc_comp=[];
        ebc_val=ebc_fenids*0;
        u   = set_ebc(u, ebc_fenids, ebc_fixed, ebc_comp, ebc_val);
        ebc_fenids=fenode_select (fens,struct('box',[a/2 a/2 0 a/2 0 h],'inflate',h/100));
        ebc_fixed=ebc_fenids*0+1;
        ebc_comp=[];
        ebc_val=ebc_fenids*0;
        u   = set_ebc(u, ebc_fenids, ebc_fixed, ebc_comp, ebc_val);
        u   = apply_ebc (u);
        % Number equations
        u   = numberdofs (u);
        % Assemble the system matrix
        K = stiffness(femm, sysmat_assembler_sparse,   geom, u);
        % Load
        corn=fenode_select (fens,struct('box',[0 0 0 0 0 h],'inflate',h/100));
        fi=force_intensity(struct('magn',[0;0;-P/4/length(corn)]));
        cfes = fe_set_P1(struct('conn',corn'));
        cfemm = femm_deformation_linear(struct('material',mater,'fes',cfes,...
            'integration_rule',point_rule));
        F = distrib_loads(cfemm, sysvec_assembler, geom, u, fi, 0);

        % Solve
        u = scatter_sysvec(u, K\F);
        ucorn= gather_values(u, corn);
        nd=abs(mean(ucorn(:,3))/w_max)
        %toc
        normalized_deflections = [normalized_deflections,nd];
        neqs = [neqs u.nfreedofs];
        if graphics
            gv=graphic_viewer;
            gv=reset (gv,[]);
            draw(fes, gv, struct ('x',geom,'u', scale*u, 'facecolor','none'));
            camset(gv, [-3.3445,-4.2190,3.1679,0.3162,0.5517, -0.3039 , 0, 0,1.0000,10.3396]);
        end
    end
end
% plot(neqs, normalized_deflections,'gd-','linewidth',3); hold on
% plot(neqs, normalized_deflections,'ro-','linewidth',3); hold on
% plot(neqs, normalized_deflections,'bx-','linewidth',3); hold on
plot(neqs, normalized_deflections,'ks-.','linewidth',3); hold on
title([ 'Span to thickness ratio =' num2str(a/h) ])