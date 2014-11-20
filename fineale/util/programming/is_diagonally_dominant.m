function  boolean = is_diagonally_dominant(C)
    % Is a matrix diagonally dominant?
    %
    %     function  boolean = is_diagonally_dominant(C)
    %
    % If the absolute value of the diagonal element iis greater than the
    % sum of the absolute values of all the other elements in a row, the
    % matrix is called diagonally dominant.
    Cdiag = spdiags(C,0);
    Crest = spdiags(zeros(size(C,1),1),0,C);
    Cdiag = spdiags(Cdiag,0,size(C,1),size(C,2)) ;
    b=sum(abs(Cdiag))>sum(abs(Crest));
    boolean = isempty(find(b==0));