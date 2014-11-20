function [fens,fes] = Q4_blockx(xs,ys,options)
    % Graded mesh  of a rectangle, Q4 finite elements.
    %
    % function [fens,fes] = Q4_blockx(xs, ys, options)
    %
    % Mesh of a 2-D block, Q4 finite elements. The nodes are located at the
    % Cartesian product of the two intervals on the input.  This allows for
    % construction of graded meshes.
    %
    % xs,ys - Locations of the individual planes of nodes.
    %
    % options - structure with fields recognized by the constructor of the
    %   fe_set_Q4 object
    % 
    % Examples:  
    %     [fens,fes] = Q4_blockx(1/125*(0:1:7).^3,4+(0:2:8), []);
    %     drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
    %  
    % See also: Q4_block, fe_set_Q4
    %

    if ~isstruct(options)
        other_dimension = options; 
        clear options;
        options.other_dimension = other_dimension;
    end
    
    nL = length(xs) - 1;
    nW = length(ys) - 1;

    nnodes = (nL+1) * (nW+1);
    ncells = nL * nW;

    % preallocate node locations
    xyz = zeros(nnodes, 2);
    k = 1;
    for j = 1:(nW+1)
        for i = 1:(nL+1)
            xyz(k,:) = [xs(i) ys(j)];
            k = k + 1;
        end
    end
    % create the nodes
    fens = fenode_set(struct('xyz', xyz));

    % preallocate connectivity matrix
    options.conn = zeros(ncells, 4);
    k = 1;
    for i = 1:nL
        for j = 1:nW
            options.conn(k,:) = node_numbers(i,j,nL,nW);
            k = k + 1;
        end
    end
    % create the cells
    fes = fe_set_Q4(options);
    
    return; % block
    
end

function nn = node_numbers(i,j,nL,nW)
    
    f = (j-1) * (nL+1) + i;
    
    nn = [f, (f+1), f+(nL+1)+1, f+(nL+1)];
    
    return;
    
end
