function retobj = combine(field1,field2)
% Make a field from a combination of two fields.
% Combine two fields into one by concatenation of the data 
% arrays of the two fields.
%
% function retobj = combine(field1,field2)
%
%   Call as:
%      newf = combine(f1,f2)
%   where
%      f1   - field with nfens parameters of dimension d1
%      f2   - field with nfens parameters of dimension d2
%      newf - field with nfens parameters of dimension d1+d2
%
    dim1=get(field1,'dim');
    dim2=get(field2,'dim');
    nfens1=get(field1,'nfens');
    nfens2=get(field2,'nfens');
    if (nfens2 ~= nfens1)
        error('Fields do not match in nfens');
    end
    retobj = field(struct ('name',[get(field1,'name') get(field2,'name')], 'dim',dim1+dim2, 'nfens',nfens1));
    retobj.values = [field1.values field2.values];
    retobj.fixed_values = [field1.fixed_values field2.fixed_values];
    retobj.is_fixed = [field1.is_fixed field2.is_fixed];
    return;
end

