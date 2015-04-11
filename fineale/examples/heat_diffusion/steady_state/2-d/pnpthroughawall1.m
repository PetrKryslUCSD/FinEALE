% Solution of the same problem as in pnpSquareInSquare but with ad hoc
% Matlab code. 
function pnpSquareInSquare1
% Coordinates of nodes
x= [0,0; 0.05,0; 0.3, 0; 0,0.1; 0.05,0.1; 0.3, 0.1];
N = 6;% total number of nodes
N_f = 2;% total number of free degrees of freedom
% Mapping from nodes to  degrees of freedom
node2dof=zeros(1,N);
node2dof([2,5]) =1:N_f;% number the free degrees of freedom
node2dof(node2dof==0) =N_f+1:N;% number the fixed degrees of freedom
% Mapping degrees of freedom to nodes
dof2node(node2dof) =(1:N)';
% Mapping of free degrees of freedom to nodes
freedof2node= dof2node(1:N_f);
% Connectivity of the left region triangles
connL = [4,1,5; 2,5,1];
% Connectivity of the right region triangles
connR = [6,5,2; 2,3,6];
syms Dz real % Thickness of the slice (it cancels out  in the end)
kappaL=0.01; %  conductivity matrix
kappaR=1.5; %  conductivity matrix
pTe=-10; pTi =+10;% Boundary conditions, degrees Celsius
% fixed temperatures at the nodes
pT =sym(zeros(N,1));
pT([1,4]) =pTe;% left face
pT([3,6]) =pTi;% right face
% gradients of the basis functions with respect to the parametric coordinates
Nder= [-1,-1;1,0;0,1];
    
Kg=sym(zeros(N_f));% allocate the global conductivity matrix
Lg=sym(zeros(N_f,1));% allocate the global heat loads vector
% Loop over the triangles in the mesh
for j=1:size(connL,1) % left domain
    J=x(connL(j,:),:)'*Nder;% compute the Jacobian matrix
    Ndersp=Nder/J;%  compute the spatial gradients of the basis functions
    Ke= (det(J)/2*Ndersp*kappaL*Ndersp'*Dz);% elementwise matrix
    edofs= node2dof(connL(j,:));%element degree of freedom array
    msk=(edofs<=N_f);% Boolean mask array: only free dofs
    % Assemble elementwise conductivity matrix to the global matrix
    Kg(edofs(msk),edofs(msk))=Kg(edofs(msk),edofs(msk))+Ke(msk,msk);
    Te=pT(connL(j,:));% Temperatures at the nodes of the element
    % Assemble elementwise load vector to the global vector
    Lg(edofs(msk))=Lg(edofs(msk))-Ke(msk,~msk)*Te(~msk);
end 
for j=1:size(connR,1) % right domain
    J=x(connR(j,:),:)'*Nder;% compute the Jacobian matrix
    Ndersp=Nder/J;%  compute the spatial gradients of the basis functions
    % Elementwise conductivity matrix: note  the material orientation
    % matrix
    Ke= (det(J)/2*Ndersp*kappaR*Ndersp'*Dz);
    edofs= node2dof(connR(j,:));%element degree of freedom array
    msk=(edofs<=N_f);% Boolean mask array: only free dofs
    Kg(edofs(msk),edofs(msk))=Kg(edofs(msk),edofs(msk))+Ke(msk,msk);
     Te=pT(connR(j,:));% Temperatures at the nodes of the element
    % Assemble elementwise load vector to the global vector
    Lg(edofs(msk))=Lg(edofs(msk))-Ke(msk,~msk)*Te(~msk);
end 

% Solve for the global temperatures at the free degrees of freedom
Tg=double (Kg\Lg)
% Initialize the array of temperatures at all the nodes from the fixed
% temperatures 
T=double (pT);
% Distribute (scatterer) the temperatures at the free degrees of freedom
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
xlabel x, ylabel y, zlabel T 
end

