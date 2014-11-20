% Static deflection of pin-pin prestressed wire with uniform load.
P= 4; % prestressing force
L =6;  % length of the span
q=0.1; % distributed load magnitude
N=3; % total number of nodes
N_f = 2; % number of free degrees of freedom
x = [0; L/2; L];% Array of node coordinates
conn= [1,2; 2,3];% Element definitions (connectivity)
% Degree-of-freedom numbers: the free first, followed by all the fixed
node2dof = [3; 1; 2];
% Mapping degrees of freedom to nodes
dof2node(node2dof) =(1:N)';
% Mapping of free degrees of freedom to nodes
freedof2node= dof2node(1:N_f);
% All deflections at all nodes. Even the degrees of freedom at supported 
% nodes  are included.
w=zeros(N,1); % note: fixed nodes have zero deflection

% Compute the global stiffness matrix of the structure
Kg=zeros(N_f,N_f);
for e=1:size(conn,1) % loop over all elements
    h =x(conn(e,2))-x(conn(e,1)); % length of the element
    Ke =P/h*[1,-1;-1,1]; % element stiffness matrix
    edofs =node2dof(conn(e,:)); % element degree-of-freedom array
    for r=1:size(Ke,1) % loop over the rows of the element matrix
        i=edofs(r);
        if (i<=N_f) % is this a free degree of freedom?
            for c=1:size(Ke,2) % loop over the columns of the element matrix
                j=edofs(c);
                if (j<=N_f) % is this a free degree of freedom?
                    Kg(i,j)=Kg(i,j)+Ke(r,c);  % assemble
                end
            end
        end
    end
end

% Compute the global load vector
Lg=zeros(N_f,1);
for e=1:size(conn,1) % loop over all elements
    h =x(conn(e,2))-x(conn(e,1)); % length of the element
    Le =q*h/2*[1;1]; % element load vector
    edofs =node2dof(conn(e,:)); % element equation array
    for r=1:size(Ke,1) % loop over the rows of the element vector
        i=edofs(r);
        if (i<=N_f) % is this a free degree of freedom?
            Lg(i)=Lg(i)+Le(r);  % assemble
        end
    end
end

% deflections for the free degrees of freedom are computed
wg=Kg\Lg;
% deflections are set at the nodes with free degrees of freedom 
w(freedof2node)=wg;

% Presentation of the deflection
figure; hold on
for e=1:size(conn,1) % loop over all elements
    line('XData',x(conn(e,:)),'YData',(2*P/q/L^2)*w(conn(e,:)),'linewidth',3);
end
grid on
set(gca,'View', [0,-90]);
xlabel('Node coordinate x','fontsize',14)
ylabel('Node deflection (2*P/q/L^2)*w','fontsize',14)
