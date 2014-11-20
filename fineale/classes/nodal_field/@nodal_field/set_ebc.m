function retobj = set_ebc (self, fenids, is_fixed, comp, val)
% Set the EBCs (essential boundary conditions).
%
% function retobj = set_ebc (self, fenids, is_fixed, comp, val)
%
% To actually apply the EBCs (that is to copy the fixed values
% into the appropriate places in the values array), use apply_ebc().
%    Call as:
%       retobj = set_ebc (self, fenids, is_fixed, comp, val)
%    where
%       self           - field
%       fenids         - array of N node identifiers
%       is_fixed  - an array of N Boolean values which say whether the component
%                        is being fixed (value ~=0), or whether
%                        it is to be free (value ==0).
%       comp           - an array of N components to which to apply the EBCs; comp
%                        may be supplied as '[]', in which case the
%                        fixed values apply to all components.
%       val            - an array of N values of fixed component magnitudes,
%                        length(val) == length(fenids)
%
% Note:  Any call to set_ebc() that would change the current assignment 
% which degrees of freedom are free and which are fixed invalidates the
% current degree-of-freedom numbering. In such a case this method sets 
%    nfreedofs = 0; and 
%    dofnums=dofnums*0+sysmat_assembler.invalid_dofnum; 
%
    previous_is_fixed=self.is_fixed;
    Maxl =max([length(fenids), length(is_fixed), length(comp), length(val)] );
    if (length(fenids)==1)&&(Maxl>1)
        fenids=repmat(fenids,Maxl,1);
    end
    if (length(is_fixed)==1)&&(Maxl>1)
        is_fixed=repmat(is_fixed,Maxl,1);
    end
    if (length(comp)==1)&&(Maxl>1)
        comp=repmat(comp,Maxl,1);
    end
    if (length(val)==1)&&(Maxl>1)
        val=repmat(val,Maxl,1);
    end
    [nfens,dim] = size(self.values);
    n = length(fenids);
    for i=1:n
        k = fenids(i);
        if (isempty(comp))
            m = 1:dim;
        else
            m = comp(i);
            if m>dim
                error(['Component: ' num2str(m) ' versus dim=' num2str(dim)]);
            end
        end
        self.is_fixed(k,m) = is_fixed(i);
        if (self.is_fixed(k,m))
            self.fixed_values(k,m) = val(i);
        end
    end
    % Any call to set_ebc() which changes which degrees of 
    %     freedom are free and which are fixed
    %     invalidates the current equation numbering.
    changed=find(self.is_fixed~=previous_is_fixed);
    if (~ isempty(changed))
        self.nfreedofs = 0;
        self.dofnums=zeros(size(self.dofnums))+sysmat_assembler_base.fixed_dofnum;
    end
    % Return the modified object
    retobj = self;
end
