function [fens,fes] = T10_quarter_cylinder_n(Radius, Length, nperradius, nL)
% Four-node tetrahedron mesh of one quarter of solid  cylinder with given number of edges per radius.
%
% Make T10 mesh of one quarter of a solid  cylinder with a given number of
% element edges per radius.
% 
% function  [fens,fes] = T10_quarter_cylinder_n(Radius, Length, nperradius, nL)
%   
% Output:
% fens= finite element node set
% fes = finite element set
%
% Examples: 
%     [fens,fes] = T10_quarter_cylinder_n(1.0, 3.0, 2, 3);
%     drawmesh({fens,fes},'fes','facecolor','red'); hold on
%
% See also: H8_cylinder
%
tol=min( [Length/nL,Radius/2/nperradius] )/100;
xyz=[0,0,0; Radius,0,0; Radius/sqrt(2),Radius/sqrt(2),0; 0,Radius,0; 0,0,Length; Radius,0,Length; Radius/sqrt(2),Radius/sqrt(2),Length; 0,Radius,Length];
[fens,fes] = H8_hexahedron(xyz,nperradius,nperradius,nL);
orientation ='b';
if (strcmp(orientation,'a'))
        t4ia = [1, 8, 5, 6; 3, 4, 2, 7; 7, 2, 6, 8; 4, 7, 8, 2; 2, 1, 6, 8; 4, 8, 1, 2];
        t4ib = [1, 8, 5, 6; 3, 4, 2, 7; 7, 2, 6, 8; 4, 7, 8, 2; 2, 1, 6, 8; 4, 8, 1, 2];
    elseif (strcmp(orientation,'b'))
        t4ia = [2,7,5,6; 1,8,5,7; 1,3,4,8; 2,1,5,7; 1,2,3,7; 3,7,8,1];
        t4ib = [2,7,5,6; 1,8,5,7; 1,3,4,8; 2,1,5,7; 1,2,3,7; 3,7,8,1];
    elseif (strcmp(orientation,'ca'))
        t4ia = [8, 4, 7, 5; 6, 7, 2, 5; 3, 4, 2, 7; 1, 2, 4, 5; 7, 4, 2, 5];
        t4ib = [7, 3, 6, 8; 5, 8, 6, 1; 2, 3, 1, 6; 4, 1, 3, 8; 6, 3, 1, 8];
    elseif (strcmp(orientation,'cb'))
        t4ia = [7, 3, 6, 8; 5, 8, 6, 1; 2, 3, 1, 6; 4, 1, 3, 8; 6, 3, 1, 8];
        t4ib = [8, 4, 7, 5; 6, 7, 2, 5; 3, 4, 2, 7; 1, 2, 4, 5; 7, 4, 2, 5];
    end
 conns=zeros(6*count(fes),4);
gc=1; ix=1;
for i=1:nperradius
    for j=1:nperradius
        for k=1:nL
            nn=fes.conn(ix,:);
            if (mod (sum( [i,j,k] ),2)==0)
                t4i =t4ib;
            else
                t4i =t4ia;
            end
            for r=1:size(t4i,1)
                conns(gc,:)=nn(t4i(r,:));
                gc=gc+1;
            end
            ix=ix+1;
        end
    end
end
fes =fe_set_T4(struct ('conn',conns(1:gc-1,:)));
[fens,fes] = T4_to_T10(fens,fes);;

bfes=mesh_boundary(fes,[]);
z1=fe_select(fens,bfes,struct('facing',true, 'direction',[0,0,-1],'tolerance',0.999));
z2=fe_select(fens,bfes,struct('facing',true, 'direction',[0,0,+1],'tolerance',0.999));
x1=fe_select(fens,bfes,struct('facing',true, 'direction',[-1,0,0],'tolerance',0.999));
y1=fe_select(fens,bfes,struct('facing',true, 'direction',[0,-1,0],'tolerance',0.999));
round1 =setdiff(1:count(bfes),[x1,y1,z1,z2] );
cn=connected_nodes(subset(bfes,round1));
for  j=1:length(cn)
    fens.xyz(cn(j),1:2)=fens.xyz(cn(j),1:2)/norm(fens.xyz(cn(j),1:2))*Radius;;
end

% drawmesh({fens,fes},'fes','facecolor','red'); hold on
end