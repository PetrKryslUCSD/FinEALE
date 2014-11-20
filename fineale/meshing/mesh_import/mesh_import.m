function [fens,fes]=mesh_import(filename,xyzscale,saveReadData)
% Import mesh. 
%
% function [fens,fes]=mesh_import(filename,xyzscale,saveReadData)
%
% General  import function to function as a front-and to import functions 
% for individual file formats.
%
% Limitations: Individual import functions may have limitations, please 
% refer to the "See also" section.
% 
% Arguments:
% filename= filename, either ABAQUS (.INP), Cosmos (.geo), Comsol (.mphtxt), 
%           or fineale  (.mat) filename.
% xyzscale= scale the coordinates by this number when they are read in, optional
% saveReadData= boolean, should the data of the saved in a format that allows 
%           for speedy reading of the same data later on? Optional.
% 
% Output:
% fens= finite element node set
% fes = finite element set
%
%
% See also: Abaqus_mesh_import, Comsol_mesh_import, geo_mesh_import, native_mesh_import, Ansys_mesh_import
%    
   if (~exist('xyzscale','var'))
        xyzscale = 1.0;
    end
    if (~exist('saveReadData','var'))
        saveReadData = ~true;
    end
    
    [pathstr, name, ext] = fileparts(filename) ;
    
    if (strcmpi(ext,'.inp'))
        [fens,fes]=Abaqus_mesh_import(filename,xyzscale,saveReadData);
    elseif (strcmpi(ext,'.geo'))
        [fens,fes]=geo_mesh_import(filename,xyzscale,saveReadData);
    elseif (strcmpi(ext,'.mphtxt'))
        [fens,fes]=Comsol_mesh_import(filename,xyzscale,saveReadData);
    elseif (strcmpi(ext,'.mat'))
        [fens,fes]=native_mesh_import(filename,xyzscale,saveReadData);
    else
        error( ['No import handler for extension ' ext]);
    end
end

