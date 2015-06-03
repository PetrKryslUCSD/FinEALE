function [fens,fes,meshinfo] = fmesher_mesh_import (meshpath)
% Import mesh created using fmesher from the xfemm package (which is a
% wrapper for Triangle)
%
% Syntax
% 
% [fens,gcells,meshinfo] = fmesher_mesh_import (filename, xyzscale, saveReadData)
%
% Inputs
%
%  meshpath - string representing the path to the mesh files saved on disk.
%    fmesher saves the mesh in several files, with the same basename. The
%    path supplied to loadmesh_mfemm is the path to any of these files but
%    without the extension. For example:
%
%    /home/myusername/somedirectoy/meshfile
%
%    where the directory /home/myusername/somedirectoy/ contains the files
%    meshfile.edge, meshfile.ele, and meshfile.node generated using
%    fmesher.
%
% 

    
    % read in the nodes
    [fid_node, delobj_node] = safefopen ([meshpath, '.node']);
    
    C = textscan (fid_node, '%d\t%d\t%d\t%d', 1, 'CommentStyle', '#');
    
    meshinfo.NNodes = C{1};
    meshinfo.NDimensions = C{2};
    meshinfo.NAttributes = C{3};
    meshinfo.HasBoundaryMarkers = C{4} == 1;
    
    if meshinfo.HasBoundaryMarkers
        formatstr = '%d\t%f\t%f\t%d';
    else
        formatstr = '%d\t%f\t%f';
    end
    
    C = textscan (fid_node, formatstr, meshinfo.NNodes, 'CollectOutput', true, 'CommentStyle', '#');
    
    fens = fenode_set (struct ('id', (1:size(C{2},1))', 'xyz', C{2}));
    
    % load the elements
    [fid_ele, delobj_ele] = safefopen ([meshpath, '.ele']);
    
    C = textscan (fid_ele, '%d\t%d\t%d', 1, 'CommentStyle', '#');
    
    meshinfo.NElements = C{1};
    meshinfo.NNodesPerElement = C{2};
    meshinfo.NElementAttributes = C{3};
    
    formatstr = ['%d\t', ...
                 repmat('%d\t', 1, meshinfo.NNodesPerElement), ...
                 repmat('%d\t', 1, meshinfo.NElementAttributes)];
              
    C = textscan (fid_ele, formatstr, meshinfo.NElements, 'CollectOutput', true, 'CommentStyle', '#');
    
    fes = fe_set_T3 ( struct ( 'conn', C{1}(:,2:(1+meshinfo.NNodesPerElement)) + 1, ...
                               'label', C{1}(:,(1+meshinfo.NNodesPerElement+1):end) ....
                             )  ...
                    );
    
end

function [fid, delobj] = safefopen (filepath)

    fid = fopen (filepath);
    
    delobj = onCleanup (@() fclose (fid));

end