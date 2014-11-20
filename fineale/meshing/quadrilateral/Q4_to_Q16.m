% Convert a mesh of Quadrilateral Q4 to Quadrilateral Q16.
%
% function [nfens,nfes] = Q4_to_Q16(fens,fes,options)
%
% Examples:
%     R=8.5;
%     [fens,fes]=Q4_sphere(R,0,1.0);
%     [fens,fes] = Q4_to_Q16(fens,fes,struct('other_dimension',1));
%     fens= onto_sphere(fens,R,[]);
%     drawmesh({fens,fes},'nodes','fes','facecolor','y', 'linewidth',2); hold on
%
% See also: fe_set_Q16

function [nfens,nfes] = Q4_to_Q16(fens,fes,options)
if (isempty(options))
	options=1.0;
end
	if ~isstruct(options)
        other_dimension = options; clear options;
        options.other_dimension = other_dimension;
    end
    if (count(fes)>100)
        box=[Inf,-Inf,Inf,-Inf];
        conns = fes.conn;
        axyz=fens.xyz;
        for i= 1:count(fes)
            conn=conns(i,:);
            for k=1:length(conn)
                xyz=axyz(conn(k),:);
                box([1,3])=min([box([1,3]);xyz]);
                box([2,4])=max([box([2,4]);xyz]);
            end
        end
        dims=(box([2,4])-box([1,3]));
        sbox=box;
        if (dims(1)>dims(2)) 
            Direction =   1;
        else
            Direction =  2;
        end
        if (Direction==1)
            sbox(2)=(box(1)+box(2))/2;
        else
            sbox(4)=(box(3)+box(4))/2;
        end
        celllist1 = fe_select(fens, fes, struct('box',sbox));
        celllist2=setdiff((1:count(fes)),celllist1);
%         disp([num2str(length(celllist1)) ',' num2str(length(celllist2))])
        if (length(celllist1)*length(celllist2)==0)
            [nfens,nfes] = do_It(fens,fes,options);
            return
        end
        [fens1,fes1] = Q4_to_Q16(fens,subset(fes,celllist1),options);
        [fens2,fes2] = Q4_to_Q16(fens,subset(fes,celllist2),options);
        [nfens,g1,g2] = merge_meshes(fens1, fes1, fens2, fes2, min(dims)/count(fes)/1000);
        nfes= cat(g1,g2);
    else
        [nfens,nfes] = do_It(fens,fes,options);
    end
end

function [fens,fes] = do_It(fens,fes,options)
    if ~isstruct(options)
        options = struct('conn',[]);
    end
    xyz1=fens.xyz;;
    allxyz= zeros(count(fes)*16,size(xyz1,2));
    nconns=zeros(count(fes),16);
    conns = fes.conn;
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
    g4=fe_set_Q4(struct('conn',(1:4)));
    nnfens=0;
    for i= 1:count(fes)
        conn = conns(i,:);
        xy4=zeros(4,size(xyz1,2));
        for k=1:length(conn)
            xy4(k,:)=xyz1(conn(k),:);
        end
        box=[min(xy4(:,1)),max(xy4(:,1)),min(xy4(:,2)),max(xy4(:,2))];
        tol=abs(min(box([2,4])-box([1,3])))/1000+eps;
        fenl=private_fenode_select_box(box, allxyz(1:nnfens,:),  tol);
        nconn =zeros(1, 16);
        for r=1:16
            N=bfun(g4,p(r,:));
            xyz=N'*xy4;
            pfenl=private_fenode_select(xyz, allxyz(1:nnfens,:),  tol);
            if (~isempty(pfenl))
                if (length(pfenl)~=1),error('more Than one node'),end
            else
                nnfens=nnfens+1;
                allxyz(nnfens,:)=xyz;
                pfenl=nnfens;
            end
            nconn(r)=pfenl;
        end
        nconns(i,:) =nconn;
    end
    fens=fenode_set(struct('xyz',allxyz(1:nnfens,:)));
    options.conn=nconns;
    fes=fe_set_Q16(options);
end

function nodelist = private_fenode_select(xyz, allxyz,  tolerance)
    xyzd=abs([allxyz(:,1)-xyz(1),allxyz(:,2)-xyz(2)]);
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
end

