function retobj = plus (A,B)
% Add a scalar to a field, or add two fields together.
%
% function retobj = plus (A,B)
%
%    Call as:
%      resfld=plus(field,double) or resfld=plus(double,field), or resfld=plus(field1,field2)
%     or
%      resfld=fld+d, resfld=d+fld, resfld=fld1+fld2
% *Both* the values attribute and the
% fixed_values attribute are modified.
%
   if     (isa(A,'nodal_field') & isa(B,'double'))
        A.values = A.values + B;
%         A.fixed_values = A.fixed_values + B;
        retobj = A;
    elseif (isa(B,'nodal_field') & isa(A,'double'))
        B.values = B.values + A;
%         B.fixed_values = B.fixed_values + A;
        retobj = B;
    elseif (isa(B,'nodal_field') & isa(A,'nodal_field'))
        B.values = B.values + A.values;
        B.fixed_values = B.fixed_values + A.fixed_values;
        retobj = B;
    else
        error('Undefined operands!');
    end
    return;
end

