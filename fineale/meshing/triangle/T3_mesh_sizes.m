function [hcurs, hests] =T3_mesh_sizes(curconn,x,targeterr,elerrs,convergence_rate)
% Estimate T3 mesh size from elementwise errors and the convergence rate.  S
% 
% function [hcurs, hests] =T3_mesh_sizes(curconn,x,targeterr,elerrs,convergence_rate)
% 
% curconn= current connectivity,
% x= current locations of nodes,
% targeterr= target error,
% elerrs= element wise errors,
% convergence_rate= estimated convergence rate

    hcurs=zeros(size(x,1),1)+Inf;
    hests=zeros(size(x,1),1)+Inf;
    for j=1:size(curconn,1)
        conn=curconn(j,:);
        V=x(conn(2),:)-x(conn(1),:);
        W=x(conn(3),:)-x(conn(1),:);
        A=(V(1)*W(2)-V(2)*W(1))/2;
        hcur=sqrt(A*4/sqrt(3));% estimate current element size
        errcur=elerrs(j);% this is the current element error
        % This is the formula from the textbook
        hest=hcur*(targeterr/errcur)^(1/convergence_rate);
        for k=1:length(conn)
            hcurs(conn(k))=min([hcurs(conn(k)),hcur]);
            hests(conn(k))=min([hests(conn(k)),hest]);
        end
    end
end