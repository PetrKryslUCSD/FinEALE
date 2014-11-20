%% Generation of the mesh of an elliptical-cross-section cylinder.
%

%% Description
%
% Determine the central transverse displacement in a simply-supported seven
% layer symmetric strip with a central line load. A 0/90/0/90/0/90/0
% material lay-up is specified with the center ply being four times as
% thick as the others.

%%
% Reference: NAFEMS Report R0031, Test No.1, 17-Dec-1998.


%%
% The plate is discretized with solid serendipity quadratic hexahedral
% elements. Because of the symmetries of the geometry and load, only the
% first-quadrant   (in XY) quarter of the plate is modeled.
%%
% The coordinate system is centered at point E (at the difference with
% respect to the original benchmark definition).  The  load is applied
% along a curve passing through point C. The simple support is applied
% along the curve passing through point B.

%%
%
% <html> <table border=0><tr><td> <img src="../docs/pub_R0031NAFEMS.jpg"
% width="50%"> </td></tr> <tr><td>Figure 1. Definition of the geometry of
% the thick elliptical plate</td></tr> </table> </html>

%% 
% We realize the simple supports along the lines  A, B and the line load at
% point C  are illegal from the point of view of convergence.  No
% convergence can be hoped for as the stress underneath the load and above
% the simple supports  is infinite in the limit (these locations are stress
% singularities).   However, for relatively coarse meshes the results away
% from the singularities are still meaningful. 
%% 
% The target quantities are displacement at the bottom surface at point E,
% the tensile axial stress at the same point,  and of the transverse shear
% stress at point D  in between the bottommost two layers (See figure 1).

%% Solution
%

function pub_R0031NAFEMS_maincode
    u= physical_units_struct;
    %%
    % The material is orthotropic, the same in all seven layers (the
    % orientation of the material is different ddepending on the layer, of
    % course).
    E1=290e3*u.MEGA*u.PA; E2=6e3*u.MEGA*u.PA; E3=E2; % Using Graphite ...
    ... Epoxy page 65 calculated G23
    G12=5e3*u.MEGA*u.PA; G13=5e3*u.MEGA*u.PA;  G23=118e3*u.MEGA*u.PA;
    nu12= 0.4; nu13= 0.02; nu23= 0.3;
    
    %%
    % The geometry of the strip.
    atilde=3*u.MM;% Length
    btilde=atilde/2;
    ell=5*atilde;
    
    %%
    % The line load is in the negative Z direction.
    q0 = -10*u.NT/u.MM;% find load
    
    %%
    % Here we define the layout and the thicknesses of the layers.
    angles =[45];
    
    
    %%
    % The mesh is created by creating a cylinder of H8 elements.
    nL=4; nperradius=2;
    %%
    % Define the geometrical tolerance using the minimal dimension in the
    % model.
    tolerance =min(atilde)/max(nL)/100;
    
    %%
    % The nodes must be located so that the simple support can be applied to an entire row of nodes.
    [fens,fes] = H8_cylinder_n(atilde, ell, nperradius, nL);
    [fens,fes] = H8_to_H64(fens,fes);
    % Right now the boundary faces are all flat.  We need to make them
    % curved. Now move all the nodes On the "cylindrical" boundary to be
    % at a given radius.
    bg=mesh_boundary(fes);
    l=[fe_select(fens,bg,struct ('facing', true, 'direction', @(x)([x(1:2),0]),'tolerance',0.01))];
    %     In the list |l| we have the faces of the "cylindrical" boundary.
    %       These are the nodes on that piece of surface.
    n=connected_nodes (subset(bg,l));
    for k=1:length(n)
        j=n(k);
        fens.xyz(j,1:2)=atilde*fens.xyz(j,1:2)/norm(fens.xyz(j,1:2));
    end
    %     Now scale the Y direction to create an elliptical cross-section.
    fens.xyz(:,2)=btilde/atilde*fens.xyz(:,2);
    %%
    % Plot the mesh.
    gv=drawmesh( {fens,fes...
        },'fes', 'facecolor','r');
    draw_axes(gv,struct('length',1.5*atilde));
    labels
    
end