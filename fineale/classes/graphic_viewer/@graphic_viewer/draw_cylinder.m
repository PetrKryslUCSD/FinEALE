function draw_cylinder(self, c1, c2, radius1, radius2, context)
% Draw cylinder.
%
% function draw_cylinder(self, c1, c2, radius1, radius2, context)
%
% required arguments:
% c1,c2= coordinates of the centers of the bases; these two points must not
%     coincide!
% radius1, radius2= radii of the bases
% context = struct with following optional fields
%    facecolor = color if polygons should be drawn in solid color,
%    tessel = number of polygons along the circumference;
%    fill=string: either 'none' or 'interp' or 'flat'
%    length_units= number, representing the length units in which the
%         graphic should be displayed
%
dim=length(c1);
if dim < 3
	c1(end+1)=0; c2(end+1)=0;
end
if dim < 2
	c1(end+1)=0; c2(end+1)=0;
end
if size(c1) == [1 3]
    c1=c1';
end
if size(c2) == [1 3]
    c2=c2';
end
if isfield(context,'tessel')
    tessel= context.tessel;
else
    tessel= 5;
end
[x,y,z]=cylinder([radius1,radius2],tessel);
z=z*norm(c2-c1);
p=surf2patch(x,y,z);
p.vertices=[p.vertices(:,1)  p.vertices(:,2)  p.vertices(:,3)];
T=transf(c1,c2);
p.vertices = p.vertices * T';
p.vertices=[p.vertices(:,1) + c1(1) ...
    p.vertices(:,2) + c1(2) ...
    p.vertices(:,3) + c1(3)];
if isfield(context,'facecolor')
    facecolor= context.facecolor;
else
    facecolor='yellow';
end
if isfield(context,'length_units')
  p.vertices=p.vertices/context.length_units;
end
h=patch('Vertices',p.vertices,'Faces',p.faces,'FaceColor',facecolor,'EdgeColor','none');
if isfield(context,'facealpha')
    set (h,'FaceAlpha', context.facealpha);
end
return;

function T=transf(c1,c2)
a=(c2-c1);
a=a/norm(a);
b=a*0; b(3)=1;
if ((norm(dot(a,b))) > 0.99999)
    b(3)=0; b(2)=1;
end
A=skewmat(a);
c=A*b;
c=c/norm(c);
b=skewmat(c)*a;
T=[b c a]; % note: cylinder() constructs it w/ axis parallel to the Z-axis!
return;


