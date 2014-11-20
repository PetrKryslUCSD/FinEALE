function [fens,fes]=native_mesh_import(filename,xyzscale,saveReadData)
% Import mesh from the native-format mesh file.
%
% function [fens,fes]=native_mesh_import(filename,xyzscale,saveReadData)
%
% Arguments:
% filename= filename,
% xyzscale= ignored,
% saveReadData= ignored
% 
% Output:
% fens= finite element node set
% fes = finite element set
%
% The native format simply stores fens,fes as the finite element node set
% and the finite elements set.
% 
% See also: mesh_import
%

% Always ignore    
    xyzscale = 1.0;
    
% Always ignore    
    saveReadData = ~true;
    
    [pathstr, name, ext] = fileparts(filename) ;
    if (isempty(pathstr)),pathstr ='.';end
    saveFile = [pathstr filesep name ,'_xyzconn' '.mat'];
    
    data=load(filename);
    
     % Create output arguments. First the nodes
    clear fens
    fens=data.fens;
    
    % Now the geometric cells for each element set
    clear fes
    fes=data.fes;
end
 
