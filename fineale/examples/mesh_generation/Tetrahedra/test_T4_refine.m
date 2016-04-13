% Test the refinement of a mesh of T4 elements.
% The initial mesh consists of elements labeled either 0 or 1; this script
% checks that the refined mesh maintains this labeling.
A=3.0;
B=1.2;
C=2;
[fens,fes] = T4_blocka(A, B, C, 4, 3, 2);
fes.label=0;% The elements are initially labeled all with 0
% We are going to select some elements and label them 1
l1=fe_select(fens,fes,struct('box',bounding_box([0,0,0; A/2, B, C]),'inflate',A/1e6));
fes.label(l1)=1;

% Refine the mesh. This should maintain the labeling of the individual connectivities.
[fens,fes] = T4_refine(fens,fes);

% Render in a graphic display.  The elements labeled 0 will be red, and the
% elements labeled 1 will be yellow.
l0=fe_select(fens,fes,struct('label',0));
gv=drawmesh({fens,subset(fes,l0)},'fes','facecolor','red')
l1=fe_select(fens,fes,struct('label',1));
gv=drawmesh({fens,subset(fes,l1)},'gv',gv,'fes','facecolor','y')

% Now concatenate the finite element set split into two pieces.  This
% checks that the labels are maintained during concatenation.
fes=cat(subset(fes,l1),subset(fes,l0));

pause(1)
%Render again in a graphic display. This time  The elements labeled 0 will be green, and the
% elements labeled 1 will be blue.
l0=fe_select(fens,fes,struct('label',0));
gv=drawmesh({fens,subset(fes,l0)},'fes','facecolor','g')
l1=fe_select(fens,fes,struct('label',1));
gv=drawmesh({fens,subset(fes,l1)},'gv',gv,'fes','facecolor','b')


