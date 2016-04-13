% Test the refinement of a mesh of T4 elements.
A=3.0;
B=1.2;
C=2;
[fens,fes] = T4_blockb(A, B, C, 4, 6, 2);
fes.label=2;% The elements are initially labeled all with 2
% We are going to select some elements and label them 1
l1=fe_select(fens,fes,struct('box',bounding_box([0,0,0; A, B/2, C]),'inflate',A/1e6));
fes.label(l1)=1;

% Create a cell array of finite element sets
fesc{1}=subset(fes,fe_select(fens,fes,struct('label',1)));
fesc{2}=subset(fes,fe_select(fens,fes,struct('label',2)));

% Create a single set that consists  of the sets stored in the cell array
fes=cat(fesc{1},fesc{2});

% Refine the mesh. This should maintain the labeling of the individual connectivities.
[fens,fes] = T4_refine(fens,fes);

% Split the resulting mesh into two sets
fesc{1}=subset(fes,fe_select(fens,fes,struct('label',1)));
fesc{2}=subset(fes,fe_select(fens,fes,struct('label',2)));

% Render in a graphic display.  The elements labeled 1 will be red, and the
% elements labeled 2 will be yellow.
gv=drawmesh({fens,fesc{1}},'fes','facecolor','red')
gv=drawmesh({fens,fesc{2}},'gv',gv,'fes','facecolor','y')
