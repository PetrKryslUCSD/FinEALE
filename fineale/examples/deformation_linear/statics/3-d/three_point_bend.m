%% Laminated Strip Under Three-Point Bending 
%

%%
% Link to the  <matlab:edit('pub_R0031NAFEMS') m-file>.
%

%% Description
%
% Determine the central transverse displacement in a simply-supported seven
% layer symmetric strip with a central line load. A 0/90/0/90/0/90/0
% material lay-up is specified with the center ply being four times as
% thick as the others.

%% 
% Reference:
% NAFEMS Report R0031, Test No.1, 17-Dec-1998.


%%
% The plate is discretized with solid quadratic hexahedral elements.
% Because of the symmetries of the geometry and load,
% only quarter of the plate is modeled.

%%
%
% <html>
% <table border=0><tr><td>
% <img src="../docs/pub_R0031NAFEMS.jpg" width="80%">
% </td></tr>
% <tr><td>Figure 1. Definition of the geometry of the thick elliptical plate</td></tr>
% </table>
% </html>

function pub_R0031NAFEMS
    u= physical_units_struct;
    E1=100e3*u.MEGA*u.PA; E2=5e3*u.MEGA*u.PA; E3=E2; 
    G12=3e3*u.MEGA*u.PA; G13=2e3*u.MEGA*u.PA;  G23=2e3*u.MEGA*u.PA; 
    nu12= 0.4; nu13= 0.02; nu23= 0.3;
    AB=30*u.MM; %  span between simple supports
    W=3*u.MM;% width of the plate
    q0 = -10*u.NT/u.MM;% find load
    
    nL=8; nO=1; nW=1;
    angles =[0];
    nts= 5*ones(length(angles),1);% number of elements per layer
    ts= [5]'*u.MM;% layer thicknesses
    TH=sum(ts);
    tolerance =min(ts)/max(nts)/10;
    
    wc_analytical=-1.06*u.MM;
    sigma11Eref=684*u.MEGA*u.PA;
    
    u_scale=20;
    
%% 
% The mesh is created using the composite-plate utility, making sure the
% nodes are placed at the location of the simple support by using the
% version |H8_composite_plate_x|.
    xs=linspace(-AB/2,AB/2,2*nL+1);
    ys=linspace(0,W/2,nW+1);
    [fens,fes] = H8_composite_plate_x(xs,ys,ts,nts);;
    [fens,fes] = H8_to_H20(fens,fes);
    
    region1list=[fe_select(fens,fes,struct('label', 1))];
    gv=drawmesh( {fens,subset(fes,region1list)...
           },'fes', 'facecolor','r');
    
    
    % Compose the model data
    clear model_data
    model_data.fens =fens;
    
    clear region
    region.property = 'orthotropic';
    region.E1 =E1;
    region.E2 =E2;
    region.E3 =E3;
    region.G12=G12;
    region.G13=G13;
    region.G23=G23;
    region.nu12=nu12;
    region.nu13=nu13;
    region.nu23=nu23;
    region.fes= subset(fes,region1list);
    region.integration_rule = gauss_rule (struct('dim', 3, 'order', 3));
    region.Rm =rotmat(angles(1)/180*pi* [0,0,1]);
    model_data.region{1} =region;
    
    clear essential
    essential.component= [1];
    essential.fixed_value= 0;
    essential.node_list = fenode_select(fens, ...
        struct('box', [0,0,-Inf,Inf,-Inf,Inf],  'inflate',tolerance));
    model_data.boundary_conditions.essential{1} = essential;
    
%% 
% 
    clear essential
    essential.component= [2];
    essential.fixed_value= 0;
    essential.node_list = fenode_select(fens, ...
        struct('box', [-Inf,Inf,0,0,-Inf,Inf],  'inflate',tolerance));
    model_data.boundary_conditions.essential{2} = essential;
    
    clear essential
    essential.component= [3];
    essential.fixed_value= 0;
    essential.node_list = [fenode_select(fens, ...
        struct('box', [AB/2,AB/2,-Inf,Inf,0,0],  'inflate',tolerance)),...
        fenode_select(fens, ...
        struct('box', [-AB/2,-AB/2,-Inf,Inf,0,0],  'inflate',tolerance))];
    model_data.boundary_conditions.essential{3} = essential;
    
    clear traction
    bdry_fes = mesh_boundary(fes, []);
    bcl = fe_select(fens, bdry_fes, ...
        struct ('box',[0,Inf,-Inf,Inf,TH,TH],'inflate',tolerance));
    line_fes = mesh_boundary(subset(bdry_fes,bcl), struct('other_dimension',1));
    lcl = fe_select(fens, line_fes, ...
        struct ('box',[0,0,-Inf,Inf,TH,TH],'inflate',tolerance));
    traction.fes =subset(line_fes,lcl);
    traction.traction= [0; 0; q0/2];
    traction.integration_rule =gauss_rule (struct('dim', 1, 'order', 3));
    model_data.boundary_conditions.traction{1} = traction;
    
    % Solve
    model_data =deformation_linear_statics(model_data);
    
    model_data.u_scale= u_scale;
    model_data=deformation_plot_deformation(model_data);
    
    nE=[fenode_select(fens, struct('box', [0,0,0,0,0,0],...
        'inflate',tolerance))];
    nC=[fenode_select(fens, struct('box', [0,0,0,0,TH,TH],...
        'inflate',tolerance))];
    nD=[fenode_select(fens, struct('box', [0,0,0,0,ts(1),ts(1)],...
        'inflate',tolerance))];
    Center_displacement=gather_values(model_data.u,nE);
    wc =mean(Center_displacement (:,3) )
        disp(['Center deflection=' num2str(wc/u.MM) ' mm, wc/wc_true=' num2str(wc/wc_analytical*100) '%'])
    
   
    %%
    % We are going to plot the stress using a nodal stress field.  It is
    % extracted from the quadrature points.  
    model_data.u_scale=u_scale;
    model_data.colormap=cadcolors2;
    model_data.stress_component=1;
    sigma11C=[]; sigma11D1=[]; sigma11D2=[]; sigma11E=[]; 
    function observer11(iregion, stressf,model_data)
        if (iregion==1),
            sigma11C= gather_values(stressf,nC);
            sigma11D1= gather_values(stressf,nD);
            sigma11E= gather_values(stressf,nE);
        else,
            sigma11D2= gather_values(stressf,nD);
        end
    end
    model_data.observer =@ observer11;
    model_data.use_spr=true;
    model_data.outputRm=eye(3);
    model_data=deformation_plot_stress(model_data);
     snapnow; % capture the  current image
   detail_camera =[-0.1106   -0.1459    0.1040    0.0047    0.0043   -0.0054   0  0    1.0    1.7035];;
   camget(model_data.postprocessing.gv)
   camset(model_data.postprocessing.gv,detail_camera);
     snapnow; % capture the  current image
   
    disp(['Point C sigma11=' num2str(sigma11C/(u.MEGA*u.PA)) ' MPa'])
    disp(['Point E sigma11=' num2str(sigma11E/(u.MEGA*u.PA)) ' MPa'])
    disp(['     to be compared with reference at E, sigma11=' num2str(sigma11Eref/(u.MEGA*u.PA)) ' MPa'])
    
    
     model_data.stress_component=5;
     %     sigma13C=[]; sigma13D1=[]; sigma13D2=[]; sigma13E=[];
     %     function observer13(iregion, stressf,model_data)
     %         if (iregion==1),
     %             sigma13C= gather_values(stressf,nC);
     %             sigma13D1= gather_values(stressf,nD);
     %             sigma13E= gather_values(stressf,nE);
     %         else,
     %             sigma13D2= gather_values(stressf,nD);
     %         end
     %     end
     %     model_data.observer =@ observer13;
        model_data.use_spr=true;
    model_data.outputRm=eye(3);
    model_data=deformation_plot_stress(model_data);
    snapnow; % capture the  current image
   camset(model_data.postprocessing.gv,detail_camera);
     snapnow; % capture the  current image
   
end
