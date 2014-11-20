% Direct time integration of pin-pin prestressed wire.  Trapezoidal rule.
function wireFEM3
L=7;% length of the wire
P=0.9;% prestressing force
mu= 80;% mass density per unit length
N=59;% number of nodes
%     The wire is simply supported: the first and the last node
N_f = N-2; % number of free degrees of freedom
x = linspace(0, L, N);% Array of node coordinates
conn= [(1:(N-1))',(2:(N))'];% Element definitions (connectivity)
% Degree-of-freedom numbers: the free dofs first, followed by all the fixed
node2dof = [N-1,(1:(N-2)),N]';% Remember: first and last node are pinned
% Mapping degrees of freedom to nodes
dof2node(node2dof) =(1:N)';
% Mapping of free degrees of freedom to nodes
freedof2node= dof2node(1:N_f);
    
% These arrays are indexed by the degree of freedom number.
w0 = zeros(N_f,1);% Initial condition: zero displacement
v0 = w0; v0(node2dof(floor(N/2)+1)) = 1;% Initial spike in velocity at the middle

% Compute the global stiffness and mass matrix of the structure
Kg=zeros(N_f,N_f);
Mg=zeros(N_f,N_f);
for e=1:size(conn,1) % loop over all elements
    h =x(conn(e,2))-x(conn(e,1)); % length of the element
    Ke =P/h*[1,-1;-1,1]; % element stiffness matrix
    Me =mu*h/2*[1,0;0,1]; % lumped element stiffness matrix
    edofs =node2dof(conn(e,:)); % element degree-of-freedom array
    for r=1:size(Ke,1) % loop over the rows of the element matrix
        i=edofs(r);
        if (i<=N_f) % is this a free degree of freedom?
            for c=1:size(Ke,2) % loop over the columns of the element matrix
                j=edofs(c);
                if (j<=N_f) % is this a free degree of freedom?
                    Kg(i,j)=Kg(i,j)+Ke(r,c);  % assemble stiffness
                    Mg(i,j)=Mg(i,j)+Me(r,c);  % assemble mass
                end
            end
        end
    end
end
% This is the right-hand side function for the Matlab integrator
function rhs =f(t,y)
    W=y(1:N_f); V=y(N_f+1:end);%extract from the solution vector
    rhs= [V;-Kg*W];
end
 
% The trapezoidal-rule integrator
 function [ts,ys] = trapezoidal(nsteps,tspan,w0,v0)
    ts= zeros(nsteps, 1);% preallocate the time
    ys = zeros(2*N_f,nsteps);% preallocate the output
    dt = (tspan(2)-tspan(1))/nsteps; %compute the time step
    t =tspan(1);% initialize the time
    for i=1:nsteps
        ts(i) =t;% save the time 
        ys(1:N_f,i) =w0; ys(N_f+1:end,i) =v0;% save the solution
        % update the velocity
        v1=(Mg+((dt/2)^2)*Kg)\((Mg-((dt/2)^2)*Kg)*v0-dt*Kg*w0);
        w1=w0+dt/2*(v0+v1);% update the displacement
        w0=w1;v0=v1;% store the solution for the next step
        t=t+dt;% update the time
    end
    ys=ys';% transpose as that format is expected from the output
end
%   Now call the integrator 
[ts,ys]=trapezoidal(3000,[0, 150],w0,v0);
% Plotting
function plot_history_surface(ts, ys)
    nsteps=length(ts);
    W=zeros(nsteps,N);
    W(:,freedof2node)= ys(:,1:N_f);
    surf(ts,x,W', 'FaceColor', 'interp', 'EdgeColor', 'none');
    hold on;
    camlight headlight ;
    contour3(ts,x,W', 10, 'k-');
    xlabel(' Time [seconds]')
    ylabel('Distance along the wire [length]')
    zlabel('Deflection [length]')
end
function plot_history_animation(ts,ys)
    nsteps=length(ts);
    vmin=min(min(ys(:,1:N_f)));
    vmax=max(max(ys(:,1:N_f)));
    for j=1:nsteps
        W=zeros(N,1);
        W(freedof2node)= ys(j,1:N_f)';
        plot(x,W);
        axis([0 L vmin vmax]);
        xlabel('Distance along the wire [length]')
        ylabel('Deflection [length]')
        left_handed_axes
        pause(1/N_f^2);
    end
end
function ke=history_kinetic_energy(ts,ys)
    nsteps=length(ts);
    nt=round(size(ys,2)/2);
    ke= 0*ts;
    for j=1:nsteps
        V=ys(j,nt+1:end);
        ke(j) =V*Mg*V'/2;
    end
end
function pe=history_potential_energy(ts,ys)
    nsteps=length(ts);
    nt=round(size(ys,2)/2);
    pe= 0*ts;
    for j=1:nsteps
        W=ys(j,1:nt);
        pe(j) =W*Kg*W'/2;
    end
end
function plot_energy(ts,ys)
    ke=history_kinetic_energy(ts,ys);
    pe=history_potential_energy(ts,ys);
    set(gca,'FontSize', 14)
    plot(ts,pe+ke,'b-','linewidth', 3)
    xlabel(' Time [seconds]')
    ylabel(' Total energy')
end

%     plot_history_surface(ts,ys)
%         plot_history_animation(ts,ys)
plot_energy(ts,ys)
end