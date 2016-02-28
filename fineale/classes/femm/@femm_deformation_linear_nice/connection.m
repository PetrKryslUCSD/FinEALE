% Compute the connection matrices of the individual gcells.
%
% function ems = connection(self)
%
% Return an array of the element matrices so they may be assembled.
%
function ems = connection(self)
    gcells = get(self,'gcells');
    % Get the inverse map from finite element nodes to geometric cells
    gmap=get(self.f2gmap,'gcell_map');
    % Prepare some data: connectivity and the element matrices
    conns = get(gcells, 'conn');
    Ke = cell(length(gmap),1);
    eqnums = cell(length(gmap),1);
    iem=1;
    for nix=1:length(gmap)
        gl=gmap{nix};
        if ~isempty(gl) % This node has a element patch in this block
            kconns = conns(gl,:); % connectivity of all geometric cells of the patch
            patchconn  = patch_conn (self,kconns,nix);
            eqnums{iem} = patchconn';
            Ke{iem} = ones(length(patchconn));
            iem=iem+1;
        end
    end
    ems  = elematset(struct('mat',{Ke}, 'eqnums',{eqnums}));
end
