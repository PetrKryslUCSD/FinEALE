function retobj = mtimes (A,B)
% Multiply a field by a scalar.
%
% function retobj = mtimes (A,B)
%
%    Call as:
%       resfld=mtimes(field,double) or resfld=mtimes(double,field)
% or
%       resfld=fld*d,  or resfld=d*fld
%
% To divide, multiply with (1/d).
%
% *Both* the values attribute and the
% fixed_values attribute are modified.
%
    if     (isa(A,'nodal_field') & isa(B,'double'))
        A.values = B * A.values;
        A.fixed_values = B * A.fixed_values;
        retobj = A;
    elseif (isa(B,'nodal_field') & isa(A,'double'))
        B.values = A * B.values;
        B.fixed_values = A * B.fixed_values;
        retobj = B;
    else
        error('Undefined operands!');
    end
    return;
end
