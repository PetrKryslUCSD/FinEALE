function self = numberdofs (self,varargin)
% Number the degrees of freedom.
%
% function retobj = numbereqns (self,varargin)
%
% Call as for field f:
%      f=numbereqns(f);
% or
%      f=numbereqns(f,options);
% where options is a struct with either 
% (1) nfreedofs =a scalar N which is the total number of degrees of freedom; the 
%          numbering of the remaining degrees of freedom will start at N+1;
% (2) dofnums =an array of equation numbers with size(n)==[nfens,dim];
% (3) dofnums_perm =permutation vector of length self.nfreedofs to be applied 
%          to the current numbering;
% (4) node_perm =permutation vector of length self.nfens to be applied 
%          to the current numbering;
% (5) linkages =cell array of linkages of nodes; each cell holds a structure
%          LINK such that LINK.slave is a list of nodes whose LINK.component
%          is constrained to be the same as the  LINK.component of the master 
%          node LINK.master;
% 
% Note that we're assigning to f: f changes inside this function!
%
% The free components in the field are numbered consecutively.
% No effort is made to optimize the numbering in any way. If you'd like to 
% optimize the numbering of the degrees of freedom, use the above form that sets
% the permutation of the degrees of freedom, or the permutation of the nodes.
%
    fixed_dofnum = sysmat_assembler_base.fixed_dofnum;
    [nfens,dim] = size(self.values);
    options= [];
    if nargin > 1
        options=varargin{1};
    end
    node_order=1:nfens;
    links=[];
    if ~isempty(options)
        if isfield(options, 'nfreedofs')
            % the current number of degrees of freedom is being modified 
            % for the subsequent numbering
            self.nfreedofs = options.nfreedofs;
        elseif isfield(options, 'dofnums')
            %the equation numbers are being supplied in an array
            if (size(self.dofnums)~=options.dofnums)
                error ('Mismatched array of equation numbers');
            end
            self.dofnums=options.dofnums;
            self.nfreedofs =max(max(self.dofnums));
            return;
        elseif isfield(options, 'dofnums_perm')
            %a permutation vector is supplied
            permutation=options.eqns_perm;
            k=1;
            for i=1:nfens
                for j=1:dim
                    if (~self.is_fixed(i,j))
                        self.dofnums(i,j) = permutation(k);
                        k=k+1;
                    end
                end
            end
            return
        elseif isfield(options, 'node_perm')
            %a permutation vector is supplied
            permutation=options.node_perm;
            node_order =node_order(permutation);
        elseif isfield(options, 'linkages')
            %cell array defining linkages between nodes is supplied
            links=zeros(size(self.dofnums));
            for i=1:length(options.linkages)
                link=options.linkages{i};
                if (length(link.slave)==1)&&(length(link.component)>1)
                    link.slave=repmat(link.slave,size( link.component ));
                end
                for j=1:length( link.component )
                    for k=1:length(link.slave)
                        M=min([link.slave(k), link.master]);
                        S=max([link.slave(k), link.master]);
                        links(S,link.component(j))=M;
                    end
                end
            end
        else
            error('Argument purpose not recognized');
        end
    else
        self.nfreedofs = 0;
    end
    % regular numbering
    for i=node_order
        for j=1:dim
            if (~self.is_fixed(i,j))% free degree of freedom
                if (~ isempty(links))
                    if (links(i,j)~=0)% is this a slave?
                        self.dofnums(i,j) = self.dofnums(links(i,j),j);
                    else
                        self.nfreedofs = self.nfreedofs + 1;
                        self.dofnums(i,j) = self.nfreedofs;
                    end
                else
                    self.nfreedofs = self.nfreedofs + 1;
                    self.dofnums(i,j) = self.nfreedofs;
                end
            else% fixed degree of freedom: no equation
                self.dofnums(i,j) = fixed_dofnum;
            end
        end
    end
end


