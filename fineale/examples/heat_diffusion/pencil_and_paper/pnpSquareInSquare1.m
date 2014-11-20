% Solution of the same problem as in pnpSquareInSquare but with ad hoc
% Matlab code. 
function pnpSquareInSquare1
% Coordinates of nodes
x= [-48,48; 0,48; 48,48; 0,31; -48,0;...
-31,0; 31,0; 48,0; 0,-31; -48,-48; ...
0,-48; 48,-48];
N = 12;% total number of nodes
N_f = 6;% total number of free degrees of freedom
% Mapping from nodes to  degrees of freedom
node2dof=zeros(1,N);
node2dof([9,7,4,6,5,8]) =1:N_f;% number the free degrees of freedom
node2dof(node2dof==0) =N_f+1:N;% number the fixed degrees of freedom
% Mapping degrees of freedom to nodes
dof2node(node2dof) =(1:N)';
% Mapping of free degrees of freedom to nodes
freedof2node= dof2node(1:N_f);
% Connectivity of the inner region triangles
conninner = [9,7,6; 7,4,6];
% Connectivity of the outer region triangles
connouter = [1,4,2; 4,3,2; 4,7,3; 1,6,4; 1,5,6; 5,10,6; 10,9,6;...
 10,11,9; 11,12,9; 12,7,9; 12,8,7; 8,3,7];
syms Dz real % Thickness of the slice (it cancels out  in the end)
kappainner=[2.25 0; 0 0.06]; % orthotropic conductivity matrix
kappaouter=[0.25 0; 0 0.25]; % isotropic conductivity matrix
alpha =-45;% local material orientation angle
ca=cos(2*pi/360*alpha); sa=sin(2*pi/360*alpha);
Rminner = [ca, -sa; sa, ca];% local material directions
Tbot=20; Ttop =57;% Boundary conditions, degrees Celsius
% fixed temperatures at the nodes
pT =sym(zeros(N,1));
pT([10, 11, 12]) =Tbot;% bottom horizontal face
pT([1, 2, 3]) =Ttop;% top horizontal face
% gradients of the basis functions with respect to the parametric coordinates
Nder= [-1,-1;1,0;0,1];
    
Kg=sym(zeros(N_f));% allocate the global conductivity matrix
Lg=sym(zeros(N_f,1));% allocate the global heat loads vector
% Loop over the triangles in the mesh
for j=1:size(connouter,1) % outer domain
    J=x(connouter(j,:),:)'*Nder;% compute the Jacobian matrix
    Ndersp=Nder/J;%  compute the spatial gradients of the basis functions
    Ke= (det(J)/2*Ndersp*kappaouter*Ndersp'*Dz);% elementwise matrix
    edofs= node2dof(connouter(j,:));%element degree of freedom array
    msk=(edofs<=N_f);% Boolean mask array: only free dofs
    % Assemble elementwise conductivity matrix to the global matrix
    Kg(edofs(msk),edofs(msk))=Kg(edofs(msk),edofs(msk))+Ke(msk,msk);
    Te=pT(connouter(j,:));% Temperatures at the nodes of the element
    % Assemble elementwise load vector to the global vector
    Lg(edofs(msk))=Lg(edofs(msk))-Ke(msk,~msk)*Te(~msk);
end 
for j=1:size(conninner,1) % inner domain
    J=x(conninner(j,:),:)'*Nder;% compute the Jacobian matrix
    Ndersp=Nder/J;%  compute the spatial gradients of the basis functions
    % Elementwise conductivity matrix: note  the material orientation
    % matrix
    Ke= (det(J)/2*Ndersp*Rminner*kappainner*Rminner'*Ndersp'*Dz);
    edofs= node2dof(conninner(j,:));%element degree of freedom array
    msk=(edofs<=N_f);% Boolean mask array: only free dofs
    Kg(edofs(msk),edofs(msk))=Kg(edofs(msk),edofs(msk))+Ke(msk,msk);
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
for j=1:size(connouter,1)
    patch('xdata',x(connouter(j,:),1),'ydata',x(connouter(j,:),2),...
        'zdata',T(connouter(j,:)), 'facecolor','y');;
end 
for j=1:size(conninner,1)
    patch('xdata',x(conninner(j,:),1),'ydata',x(conninner(j,:),2),...
        'zdata',T(conninner(j,:)), 'facecolor','y');;
end 
axis equal
xlabel x, ylabel y, zlabel T 
