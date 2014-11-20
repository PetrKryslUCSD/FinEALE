function val = times(A,B)
% Term by term multiplication of two fields (scalar, or dot, product).
%
% function val = times(A,B)
%
% The resulting field has as degrees of freedom the dot products of 
% the degrees of freedom arrays for each node.
    if (isa(B,'nodal_field') & isa(A,'nodal_field'))
        if (size(A.values) ~= size(B.values))
            error('Size mismatch!');
        end
        [nfens,dim] = size(B.values);
        d = 0;
        for i=1:nfens
            d = d + dot(A.values(i,:),B.values(i,:));
        end
        val = d;
    else
        error ('Invalid operands!');
    end
end
