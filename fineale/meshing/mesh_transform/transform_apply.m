function fens = transform_apply(fens, function_handle, function_data)
% Apply geometry transformation to the locations of the nodes.
%
% function fens = transform_apply(fens, function_handle, function_data)
%
% Apply geometry transformation to the locations of the nodes. 
% The function handle function_handle is of type
%    function xyz=fun(xyz, function_data)
% Where 
% xyz = location of a given node,
% function_data = anything that the function might need to do its job.
    xyz=fens.xyz;
    x=function_handle (xyz(1,:), function_data);
    xyz1 =zeros(size(xyz,1),size(x,2));
    for m=1:size(xyz,1)
        xyz1(m,:) =  function_handle (xyz(m,:), function_data);
    end
    fens.xyz=xyz1;
end
