% Plot the data as a stream tube.
%
% function streamplot(t, vdata)
%
function streamplot(t, vdata)
n=length(t);
nc=n-1;
r=0.016;
tessel=8;
mint=min(t);
maxt=max(t);
dcmap=data_colormap(mint, maxt,'colormap',jet);
for i=1:nc
    [x,y,z]=cylinder([r r],tessel);
    c1=vdata(:,i);
    c2=vdata(:,i+1);
    z=z*norm(c2-c1);
    p=surf2patch(x,y,z);
    p.vertices=[p.vertices(:,1)  p.vertices(:,2)  p.vertices(:,3)];
    T=transf(c1,c2); 
    p.vertices = p.vertices * T';
    p.vertices=[p.vertices(:,1) + c1(1) ... 
            p.vertices(:,2) + c1(2) ...
            p.vertices(:,3) + c1(3)];
    color=map_data(dcmap, t(i));
    patch('Vertices',p.vertices,'Faces',p.faces,'FaceColor',color,'EdgeColor','none');
end
return;

function T=transf(c1,c2)
if (size(c1)==[3 1])
    a=(c2-c1);
else
    a=(c2-c1)';
end
a=a/norm(a);
b=a*0; b(3)=1;
if ((norm(a-b)) < 1e-6)
    b(3)=0; b(2)=1;
end
A=skewmat(a);
c=A*b;
c=c/norm(c);
b=skewmat(c)*a;
T=[b c a]; % note: cylinder() constructs it w/ axis parallel to the Z-axis!
return;
