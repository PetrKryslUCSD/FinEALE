function model_data =renumber_mesh(model_data, method)
% Algorithm for computing the "optimal" numbering of nodes.
%
% function model_data =renumber_mesh(model_data, method)
%
% Compute an optimal numbering of the finite element nodes using an
% appropriate method. The numbering may  be then passed to
% renumbering_mesh_update() to renumber the mesh.
% S = connection matrix
% method =  either 'symrcm' or 'symamd' (optional argument, default
%           'symrcm')
%
% Output
% model_data.numbering = struct with fields P (permutation), and IP (inverse permutation)

    
    if (~exist('method','var'))
        method ='symrcm';%  default
    end
    
    % Construct the geometry field
    geom = nodal_field(struct('name',['geom'], 'dim', model_data.fens.dim, 'fens', model_data.fens));
    
    % Construct the temperature field
    temp = nodal_field(struct('name',['temp'], 'dim', 1, 'nfens', geom.nfens));
    temp = numberdofs (temp);
    
    S=  sparse(temp.nfens,temp.nfens);
    for i=1:length(model_data.region)
        region =model_data.region{i};
        if (isfield(region,'fes'))
            femm = femm_base (struct ('material',[],'fes',region.fes));
        else
            femm = region.femm;
        end
        S = S + connection_matrix(femm, sysmat_assembler_sparse, temp);
        clear region femm
    end
    
    if (strcmp(method,'symrcm'))
        % Reverse Cuthill-McKee
        p=symrcm(S);
    elseif (strcmp(method,'symamd'))
        p=symamd(S);
    else
        eval(['p=' method '(S);']);
    end
    
    clear S
    
    ip=(1:length(p));
    ip(p)=(1:length(p));
    model_data.node_perm =p; % this is the node permutation
end
