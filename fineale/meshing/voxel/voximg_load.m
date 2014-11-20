% Load a Matlab file storing a voxel image.
%
% function [img,dx,dy,dz] = voximg_load(name)
%
% The file is expected to store variables
% img = 3-D image (array), the voxel values  are arbitrary
% dx, dy, dz =dimensions of voxels 
% 
function [img,dx,dy,dz] = voximg_load(name)
     img_data= load (name);
     img=img_data.img;
     if (isfield(img_data, 'dx'))
         dx=img_data.dx;
         dy=img_data.dy;
         dz=img_data.dz;
     elseif (isfield(img_data, 'pixdim'))
         dx=img_data.pixdim(1);
         dy=img_data.pixdim(2);
         dz=img_data.pixdim(3);
     else
         dx=[];
         dy=[];
         dz=[];
     end
end