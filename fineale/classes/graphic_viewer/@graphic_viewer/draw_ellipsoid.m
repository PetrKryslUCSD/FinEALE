function draw_ellipsoid(self, c, princ_dir, radii, context)
  % Draw ellipsoid.
  %
  % function draw_ellipsoid(self, c, princ_dir, radii, context)
  %
  % required arguments:
  % c= coordinates of the center; tthis function expects to draw a single
  %      ellipsoid
  % princ_dir= matrix of principal directions in the columns;
  % radii= radii along the principle directions;
  % context = struct with following optional fields
  %    facecolor = color if polygons should be drawn in solid color,
  %    tessel = number of polygons along the circumference;
  %    length_units= number, representing the length units in which the
%         graphic should be displayed
%
  
  dim=length(c);
  if dim<3
    c(end+1)=0;
  end
  if dim<2
    c(end+1)=0;
  end
  c=reshape(c,3,1);
  if isfield(context,'tessel')
    tessel= context.tessel;
  else
    tessel= 6;
  end
  if isfield(context,'facecolor')
    facecolor= context.facecolor;
  else
    facecolor='yellow';
  end
  [x,y,z]=myellipsoid(radii(1),radii(2),radii(3),tessel);
  p=surf2patch(x,y,z);
  p.vertices=[p.vertices(:,1)  p.vertices(:,2)  p.vertices(:,3)];
  p.vertices = p.vertices * princ_dir';
  p.vertices=[p.vertices(:,1) + c(1) ...
    p.vertices(:,2) + c(2) ...
    p.vertices(:,3) + c(3)];
  if isfield(context,'length_units')
        p.vertices=p.vertices/context.length_units;
  end
  p=patch('Vertices',p.vertices,'Faces',p.faces,...
    'FaceColor',facecolor,'EdgeColor','none');
  if isfield(context,'facealpha')
    set (p,'FaceAlpha', context.facealpha);
  end
  return;

function [xx,yy,zz] = myellipsoid(r1,r2,r3,n)
% -pi <= theta <= pi is a row vector.
% -pi/2 <= phi <= pi/2 is a column vector.
theta = (-n:2:n)/n*pi;
phi = (-n:2:n)'/n*pi/2*0.999;
cosphi = cos(phi); %cosphi(1) = 0; cosphi(n+1) = 0;
sintheta = sin(theta); sintheta(1) = 0; sintheta(n+1) = 0;
xx = r1*cosphi*cos(theta);
yy = r2*cosphi*sintheta;
zz = r3*sin(phi)*ones(1,n+1);
return;
