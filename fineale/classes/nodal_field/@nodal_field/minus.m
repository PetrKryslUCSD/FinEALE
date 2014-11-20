function self = minus (A,B)
% Subtract a scalar from a field, or subtract a field from another field.
%
% function self = minus (A,B)
%
%    Call as:
%      resfld=minus(field,double) or resfld=minus(double,field), or resfld=minus(field1,field2)
%     or
%      resfld=fld-d, resfld=d-fld, resfld=fld1-fld2
%
% *Both* the values attribute and the
% fixed_values attribute are modified.
%
    if     (isa(A,'nodal_field') & isa(B,'double'))
        A.values = A.values - B;
        A.fixed_values = A.fixed_values - B;
        self = A;
    elseif (isa(B,'nodal_field') & isa(A,'double'))
        B.values = B.values - A;
        B.fixed_values = B.fixed_values - A;
        self = B;
    elseif (isa(B,'nodal_field') & isa(A,'nodal_field'))
        A.values = A.values - B.values;
        A.fixed_values = A.fixed_values - B.fixed_values;
        self = A;
    else
        error('Undefined operands!');
    end
end
