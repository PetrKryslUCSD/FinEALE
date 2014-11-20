function S = connection_matrix (self, assembler, u)
% Compute the connection matrix.
% Compute the connection matrix by computing and assembling the
% matrices of the individual FEs.
%
% function S = connection_matrix (self, assembler, u)
%
% Arguments:
%     assembler = descendent of the sysmat_assembler class
%     u=displacement field
%
% Return a matrix representing the connections between the finite element nodes
% expressed by the  finite elements.
    fes = self.fes;% grab the finite elements to work on
    % Retrieve data for efficiency
    tconns =transpose(fes.conn);% connectivity, Transposed for efficiency
    % Prepare assembler
    Kedim =fes.nfens; Ke = ones(Kedim,Kedim);
    start_assembly(assembler, Kedim, Kedim, size(tconns,2), u.nfens, u.nfens);
    % Now loop over all fes in the block
    for i=1:size(tconns,2)
        assemble_symmetric(assembler, Ke, tconns(:,i));
    end
    S = make_matrix (assembler);
end
