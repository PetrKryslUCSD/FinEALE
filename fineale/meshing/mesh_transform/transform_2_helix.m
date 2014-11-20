function fens = transform_2_helix(fens,climbPerRevolution)
% Generate the mesh of an helix. 
%
% function fens = transform_2_helix(fens,climbPerRevolution)
%
% Generate the mesh of an helix. It is
% assumed that in the X-direction the coordinate represents angle, and the
% Y and Z coordinates represent the cross-section.
%
% See also: transform_apply

    data.climbPerRevolution=climbPerRevolution;
    fens = transform_apply(fens,@xf, data);
end


function xyz = xf(xyz, data)
    xbar=xyz(1);
    ybar=xyz(2);
    zbar=xyz(3);
    xyz= [ybar*cos(-xbar),...
        ybar*sin(-xbar),...
        zbar+data.climbPerRevolution*xbar/(2*pi)];
end
