% Test the refinement of a mesh of T4 elements.
% The initial mesh consists of elements labeled either 0 or 1; this script
% checks that the refined mesh maintains this labeling.
A=3.0;
B=1.2;
C=2;
[fens,fes] = T4_blocka(A, B, C, 4, 3, 2);
fes.label=0;
l1=fe_select(fens,fes,struct('box',bounding_box([0,0,0; A/2, B, C]),'inflate',A/1e6));
fes.label(l1)=1;
[fens,fes] = T4_refine(fens,fes);

l0=fe_select(fens,fes,struct('label',0));
gv=drawmesh({fens,subset(fes,l0)},'fes','nodes','facecolor','red')
l1=fe_select(fens,fes,struct('label',1));
gv=drawmesh({fens,subset(fes,l1)},'gv',gv,'fes','nodes','facecolor','y')
