function fens = block_function_shift(fens, options)
% Shift nodes in a three-dimensional block using a displacement function.
%
% function fens = block_function_shift(fens, options)
%
% fens = array of finite element nodes
% options = struct with mandatory field
%      displacement = cell array of function handles of the signature
%           function u=f(xyz), where xyz is an array of the coordinates
%           (one triple per row).
% The locations of the nodes are only changed in the interior so that
% nodes on the boundary surface remain on the boundary.
    xyz =fens.xyz;
    Dimension =size(xyz,2);
    xyz0 = xyz;
    
    minxyz =[];
    maxxyz =[];
    for j=1:Dimension
        minxyz(j) =min(xyz(:,j));
        maxxyz(j) =max(xyz(:,j));
    end
    
    tol = [];
    if isfield(options,'tol')
        tol = options.tol;
    else
        tol = min([(maxxyz-minxyz)/1e7]);
    end
    
    for j=1:Dimension
        displacement_function =options.displacement{j};
        xyz(:,j)=xyz0(:,j) + displacement_function (xyz0);
    end

    %     make sure the nodes on the surfaces stay on the surfaces
    for j=1:Dimension
        box = [];
        for k=1:Dimension
            if j==k
                box = [box, minxyz(k), minxyz(k)];
            else
                box = [box, minxyz(k), maxxyz(k)];
            end
        end
        minxnl=fenode_select (fens,struct ('box',box,'inflate',tol));
        box = [];
        for k=1:Dimension
            if j==k
                box = [box, maxxyz(k), maxxyz(k)];
            else
                box = [box, minxyz(k), maxxyz(k)];
            end
        end
        maxxnl=fenode_select (fens,struct ('box',box,'inflate',tol));
        xyz(minxnl,j) =minxyz(j);
        xyz(maxxnl,j) =maxxyz(j);
    end
    
    %     now set the new coordinates
    fens.xyz=xyz;
end
  
