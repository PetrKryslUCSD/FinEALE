% FRF calculation.
function wireFEM4
L=7;% length of the wire
P=900;% prestressing force
mu= 80;% mass density per unit length
N=57;% number of nodes
%     The wire is on rollers at both ends.
N_f = N; % number of free degrees of freedom
x = linspace(0, L, N);% Array of node coordinates
conn= [(1:(N-1))',(2:(N))'];% Element definitions (connectivity)
% Degree-of-freedom numbers: the free dofs first, followed by all the fixed
node2dof = [(1:N)]';% All nodes are free
% Mapping degrees of freedom to nodes
dof2node(node2dof) =(1:N)';
% Mapping of free degrees of freedom to nodes
freedof2node= dof2node(1:N_f);
    
% Compute the global stiffness and mass matrix of the structure. 
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

Fg=zeros(N_f,1);
Fg(1) =1;

omegas = linspace(0.001,2,1000)*2*pi;

Us=zeros(size(omegas));

for i=1:length(omegas)
    omega =omegas(i);
    U = (-omega^2*Mg+Kg)\Fg;
    Us(i) = U(N_f);
end

plot(omegas,abs(Us))
set(gca,'ylim',[0, 0.45 ])
grid on
xlabel('Angular frequency')
ylabel('FRF')


% Natural frequencies
k=0:1:5;
omega_n =pi/L*k*sqrt(P/mu)
end