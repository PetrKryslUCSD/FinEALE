function [nfens,nfes] = H8_to_H64(fens,fes)
% Convert a mesh of hexahedra H8 to hexahedra H64.
%
% function [fens,fes] = H8_to_H64(fens,fes)
%
% Arguments and
% Output:
% fens= finite element node set
% fes = finite element set
%
    if (count(fes)>500)
        box=[Inf,-Inf,Inf,-Inf,Inf,-Inf];
        conns = fes.conn; 
        xyzs=fens.xyz;
        for i= 1:size(conns,1)
            conn = conns(i,:);
            for k=1:length(conn)
                xyz=xyzs(conn(k),:);
                box([1,3,5])=min([box([1,3,5]);xyz]);
                box([2,4,6])=max([box([2,4,6]);xyz]);
            end
        end
        dims=(box([2,4,6])-box([1,3,5]));
        sbox=box;
        if (dims(3)>dims(2)) && (dims(3)>dims(1))
            Direction =  3;
        elseif (dims(1)>dims(2)) && (dims(1)>dims(3))
            Direction =   1;
        else
            Direction =  2;
        end
        if (Direction==3)
            sbox(6)=(box(5)+box(6))/2;
        elseif (Direction==1)
            sbox(2)=(box(1)+box(2))/2;
        else
            sbox(4)=(box(3)+box(4))/2;
        end
        celllist1 = fe_select(fens, fes, struct('box',sbox,'anynode',1));
        celllist2=setdiff((1:count(fes)),celllist1);
        if (length(celllist1)*length(celllist2)==0)
            [nfens,nfes] = do_H8_to_H64(fens,fes);
            return
        end
        [fens1,fes1] = H8_to_H64(fens,subset(fes,celllist1));
        [fens2,fes2] = H8_to_H64(fens,subset(fes,celllist2));
        [nfens,g1,g2] = merge_meshes(fens1, fes1, fens2, fes2, min(dims)/count(fes)/1000);
        nfes= cat(g1,g2);
    else
        [nfens,nfes] = do_H8_to_H64(fens,fes);
    end
end

function [nfens,nfes] = do_H8_to_H64(fens,fes)
    %     number of vertex nodes
    p(1,:)=[-1,-1,-1];
    p(2,:)=[-0.333333333333333,-1,-1];
    p(3,:)=[0.333333333333333,-1,-1];
    p(4,:)=[1,-1,-1];
    p(5,:)=[-1,-0.333333333333333,-1];
    p(6,:)=[-0.333333333333333,-0.333333333333333,-1];
    p(7,:)=[0.333333333333333,-0.333333333333333,-1];
    p(8,:)=[1,-0.333333333333333,-1];
    p(9,:)=[-1,0.333333333333333,-1];
    p(10,:)=[-0.333333333333333,0.333333333333333,-1];
    p(11,:)=[0.333333333333333,0.333333333333333,-1];
    p(12,:)=[1,0.333333333333333,-1];
    p(13,:)=[-1,1,-1];
    p(14,:)=[-0.333333333333333,1,-1];
    p(15,:)=[0.333333333333333,1,-1];
    p(16,:)=[1,1,-1];
    p(17,:)=[-1,-1,-0.333333333333333];
    p(18,:)=[-0.333333333333333,-1,-0.333333333333333];
    p(19,:)=[0.333333333333333,-1,-0.333333333333333];
    p(20,:)=[1,-1,-0.333333333333333];
    p(21,:)=[-1,-0.333333333333333,-0.333333333333333];
    p(22,:)=[-0.333333333333333,-0.333333333333333,-0.333333333333333];
    p(23,:)=[0.333333333333333,-0.333333333333333,-0.333333333333333];
    p(24,:)=[1,-0.333333333333333,-0.333333333333333];
    p(25,:)=[-1,0.333333333333333,-0.333333333333333];
    p(26,:)=[-0.333333333333333,0.333333333333333,-0.333333333333333];
    p(27,:)=[0.333333333333333,0.333333333333333,-0.333333333333333];
    p(28,:)=[1,0.333333333333333,-0.333333333333333];
    p(29,:)=[-1,1,-0.333333333333333];
    p(30,:)=[-0.333333333333333,1,-0.333333333333333];
    p(31,:)=[0.333333333333333,1,-0.333333333333333];
    p(32,:)=[1,1,-0.333333333333333];
    p(33,:)=[-1,-1,0.333333333333333];
    p(34,:)=[-0.333333333333333,-1,0.333333333333333];
    p(35,:)=[0.333333333333333,-1,0.333333333333333];
    p(36,:)=[1,-1,0.333333333333333];
    p(37,:)=[-1,-0.333333333333333,0.333333333333333];
    p(38,:)=[-0.333333333333333,-0.333333333333333,0.333333333333333];
    p(39,:)=[0.333333333333333,-0.333333333333333,0.333333333333333];
    p(40,:)=[1,-0.333333333333333,0.333333333333333];
    p(41,:)=[-1,0.333333333333333,0.333333333333333];
    p(42,:)=[-0.333333333333333,0.333333333333333,0.333333333333333];
    p(43,:)=[0.333333333333333,0.333333333333333,0.333333333333333];
    p(44,:)=[1,0.333333333333333,0.333333333333333];
    p(45,:)=[-1,1,0.333333333333333];
    p(46,:)=[-0.333333333333333,1,0.333333333333333];
    p(47,:)=[0.333333333333333,1,0.333333333333333];
    p(48,:)=[1,1,0.333333333333333];
    p(49,:)=[-1,-1,1];
    p(50,:)=[-0.333333333333333,-1,1];
    p(51,:)=[0.333333333333333,-1,1];
    p(52,:)=[1,-1,1];
    p(53,:)=[-1,-0.333333333333333,1];
    p(54,:)=[-0.333333333333333,-0.333333333333333,1];
    p(55,:)=[0.333333333333333,-0.333333333333333,1];
    p(56,:)=[1,-0.333333333333333,1];
    p(57,:)=[-1,0.333333333333333,1];
    p(58,:)=[-0.333333333333333,0.333333333333333,1];
    p(59,:)=[0.333333333333333,0.333333333333333,1];
    p(60,:)=[1,0.333333333333333,1];
    p(61,:)=[-1,1,1];
    p(62,:)=[-0.333333333333333,1,1];
    p(63,:)=[0.333333333333333,1,1];
    p(64,:)=[1,1,1];
    g8=fe_set_H8(struct('conn',(1:8)));
    nnfens=0;
    allxyz= zeros(count(fes)*64,3);
    xyzs=fens.xyz;
    conns = fes.conn;
    nconns= zeros(count(fes),64);
    for i= 1:count(fes)
        conn = conns(i,:);
        xyz8=xyzs(conn,:);
        box=[min(xyz8(:,1)),max(xyz8(:,1)),min(xyz8(:,2)),max(xyz8(:,2)),min(xyz8(:,3)),max(xyz8(:,3))];
        tol=abs(min(box([2,4,6])-box([1,3,5])))/1000+eps;
        %         fenl=fenode_select(nfens(1:nnfens),struct('box',box,'inflate',tol));
        fenl=private_fenode_select_box(box, allxyz(1:nnfens,:),  tol);
        for r=1:64
            N=bfun(g8,p(r,:));
            xyz=N'*xyz8;
            pfenl=private_fenode_select(xyz, allxyz(1:nnfens,:),  tol);
            if (~isempty(pfenl))
                if (length(pfenl)~=1),error('more Than one node'),end
            else
                nnfens=nnfens+1;
                allxyz(nnfens,:)=xyz;
                pfenl=nnfens;
            end
            nconns(i,r)=pfenl;
        end
        %         i
    end
    nfes =fe_set_H64(struct('label',fes.label,'conn',nconns));
    allxyz=allxyz(1:nnfens,:);
    nfens=fenode_set(struct('xyz',allxyz));
end

function nodelist = private_fenode_select(xyz, allxyz,  tolerance)
    xyzd=abs([allxyz(:,1)-xyz(1),allxyz(:,2)-xyz(2),allxyz(:,3)-xyz(3)]);
    nodelist=find(sum(xyzd')'<tolerance);
end


function nodelist = private_fenode_select_box(Box, allxyz,  tolerance)
    nodelist=find(allxyz(:,1)>=Box(1));
    if (isempty(nodelist)), return, end
    nodelis1=find(allxyz(:,1)<=Box(2));
    if (isempty(nodelis1)), return, end
    nodelist =union(nodelist,nodelis1);
    if (isempty(nodelist)), return, end
    
    nodelis1=find(allxyz(:,2)>=Box(3));
    if (isempty(nodelis1)), return, end
    nodelist =union(nodelist,nodelis1);
    if (isempty(nodelist)), return, end
    nodelis1=find(allxyz(nodelist,2)<=Box(4));
    if (isempty(nodelis1)), return, end
    nodelist =union(nodelist,nodelis1);
    if (isempty(nodelist)), return, end
    
    nodelis1=find(allxyz(:,3)>=Box(5));
    if (isempty(nodelis1)), return, end
    nodelist =union(nodelist,nodelis1);
    if (isempty(nodelist)), return, end
    nodelis1=find(allxyz(nodelist,3)<=Box(6));
    if (isempty(nodelis1)), return, end
    nodelist =union(nodelist,nodelis1);
end

