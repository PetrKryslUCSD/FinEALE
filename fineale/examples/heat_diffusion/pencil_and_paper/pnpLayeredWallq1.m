% Finite Element Modeling with Abaqus and Matlab for  Thermal and 
% Stress Analysis
% (C)  2015, Petr Krysl
%
% Layered wall.  Temperature boundary condition on the left, and heat flux
% boundary condition on the right.
function pnpLayeredWallq1
% Coordinates of nodes
x= [0,0; 0.07,0; 0.3,0; 0,0.1; 0.07,0.1; 0.3,0.1];
N = 6;% total number of nodes
N_f = 4;% total number of free degrees of freedom
% Mapping from nodes to  degrees of freedom
node2dof=zeros(1,N);
node2dof([2,5,3,6]) =1:N_f;% number the free degrees of freedom
node2dof(node2dof==0) =N_f+1:N;% number the fixed degrees of freedom
% Mapping degrees of freedom to nodes
dof2node(node2dof) =(1:N)';
% Mapping of free degrees of freedom to nodes
freedof2node= dof2node(1:N_f);
% Connectivity of the left region triangles
connL = [4,1,5; 2,5,1];
% Connectivity of the right region triangles
connR = [6,5,2; 2,3,6];
%Connectivity  of the boundary L2 element  on the right
connR2=[3,6];
% syms Dz real % Thickness of the slice (it cancels out  in the end)
Dz=1.0;
kappaL=0.05; %  conductivity matrix
kappaR=1.8; %  conductivity matrix
pTe=-10; % Boundary conditions, degrees Celsius
qnbar=-30;
% fixed temperatures at the nodes
pT =sym(zeros(N,1));
pT([1,4]) =pTe;% left face
% gradients of the basis functions with respect to the parametric coordinates
gradNpar= [-1,-1;1,0;0,1];
    
Kg=sym(zeros(N_f));% allocate the global conductivity matrix
Lg=sym(zeros(N_f,1));% allocate the global heat loads vector
% Loop over the triangles in the mesh
for j=1:size(connL,1) % left domain
    J=x(connL(j,:),:)'*gradNpar;% compute the Jacobian matrix
    gradN=gradNpar/J;% compute the spatial gradients of the basis functions
    Ke= (det(J)/2*gradN*kappaL*gradN'*Dz);% elementwise matrix
    edof= node2dof(connL(j,:));%element degree of freedom array
    msk=(edof<=N_f);% Boolean mask array: only free dofs
    % Assemble elementwise conductivity matrix to the global matrix
    Kg(edof(msk),edof(msk))=Kg(edof(msk),edof(msk))+Ke(msk,msk);
    Te=pT(connL(j,:));% Temperatures at the nodes of the element
    % Assemble elementwise load vector to the global vector
    Lg(edof(msk))=Lg(edof(msk))-Ke(msk,~msk)*Te(~msk);
end 
for j=1:size(connR,1) % right domain
    J=x(connR(j,:),:)'*gradNpar;% compute the Jacobian matrix
    gradN=gradNpar/J;%  compute the spatial gradients of the basis functions
    % Elementwise conductivity matrix: note  the material orientation
    % matrix
    Ke= (det(J)/2*gradN*kappaR*gradN'*Dz);
    edof= node2dof(connR(j,:));%element degree of freedom array
    msk=(edof<=N_f);% Boolean mask array: only free dofs
    Kg(edof(msk),edof(msk))=Kg(edof(msk),edof(msk))+Ke(msk,msk);
    Te=pT(connR(j,:));% Temperatures at the nodes of the element
    if (any(~msk))
        % Assemble elementwise load vector to the global vector
        Lg(edof(msk))=Lg(edof(msk))-Ke(msk,~msk)*Te(~msk);
    end
end 

for j=1:size(connR2,1) % right bboundary
    he=norm(diff(x(connR2(j,:),:)));
    Leq=-qnbar*he*Dz/2*[1;1];
    edof= node2dof(connR2(j,:));%element degree of freedom array
    msk=(edof<=N_f);% Boolean mask array: only free dofs
    % Assemble elementwise load vector to the global vector
    Lg(edof(msk))=Lg(edof(msk))+Leq(msk);
end

% Solve for the global temperatures at the free degrees of freedom
Tg=double (Kg\Lg)
% Initialize the array of temperatures at all the nodes from the fixed
% temperatures 
T=double (pT);
% Distribute (scatter) the temperatures at the free degrees of freedom
% into the array of temperatures at the nodes 
T(freedof2node) =Tg

% Draw the temperature distribution  on the mesh
for j=1:size(connL,1)
    patch('xdata',x(connL(j,:),1),'ydata',x(connL(j,:),2),...
        'zdata',T(connL(j,:))/100, 'facecolor','y');;
end 
for j=1:size(connR,1)
    patch('xdata',x(connR(j,:),1),'ydata',x(connR(j,:),2),...
        'zdata',T(connR(j,:))/100, 'facecolor','y');;
end 
axis equal
grid on
xlabel x, ylabel y, zlabel 'T/100' 
view(3)
