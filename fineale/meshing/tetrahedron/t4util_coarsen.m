function [t,v,tmid] = t4util_coarsen(t,v,tmid,options)
% Coarsen a T4 (tetrahedral) mesh
%
% function [t,v,tmid] = t4util_coarsen(t,v,tmid,options)
%
% t = array with vertex numbers, one per tetrahedron in each row
% v= array of vertex coordinates, one per vertex in each row
% tmid = tetrahedron material identifiers, one for tetrahedron in each row
% options = structure with optional fields
%    bv=array for all vertices in the input array v.  
%           Is the vertex on the boundary? True or false. The boundary vertices 
%           are in layer 1. The vertex layer numbers increase away from 
%           the boundary,  and the bigger the vertex layer number the bigger 
%           the mesh size.
%    desired_ts=mesh size desired for vertices in layer 1, here mesh size is 
%           the length of an edge
%    stretch=the mesh size should increase by this much within one layer of  
%           elements, default is 1.25
%    nblayer=number of boundary layers where the mesh size should not increase, 
%           default is 1
%    surface_coarsening = is the coarsening intended for the interior or for 
%           the surface? default is false, which means the default 
%           is coarsening of the interior.
%    preserve_thin= when coarsening the surface, should features which are  thin 
%           be unaffected by the coarsening? Here thin means consisting of 
%           only "surface" vertices.
%    vertex_weight= weight of vertices, one per  vertex; weight <= 1.0 is ignored, 
%           but weight>1.0  is used to "stretch" edges connected  to the vertex.
%           In this way one may preserve certain edges by assigning  larger 
%           weight to their vertices. Default is vertex_weight= [] (which means 
%           ignore the vertex weights)
    desired_ts =[];
    stretch =1.25;
    bv = []; % boundary vertices are marked by true, non-boundary vertices are marked by false
    % vlayer vertex layer as number, increasing from the boundary (1)
    % The bigger the vertex layer number,  the bigger the associated mesh size
    vlayer = [];
    nblayer=1;
    surface_coarsening = false;
    preserve_thin = false;
    vertex_weight= [];
    if (exist('options','var'))
        if (isfield(options,'bv'))
            vlayer=zeros(size(options.bv,1),1,'int16');
            vlayer(find (options.bv==true))=1; % start from the boundary
        end
        if (isfield(options,'desired_ts'))
            desired_ts= options.desired_ts;
        end
        if (isfield(options,'stretch'))
            stretch= options.stretch;
        end
        if (isfield(options,'nblayer'))
            nblayer= options.nblayer;
        end
        if (isfield(options,'vertex_weight'))
            vertex_weight= options.vertex_weight;
        end
        if (isfield(options,'surface_coarsening'))
            surface_coarsening= options.surface_coarsening;
        end
        if (isfield(options,'preserve_thin'))
            preserve_thin= options.preserve_thin;
        end
    end
    
    v2t = t4util_v_to_t4_map (size(v,1),t) ;% Map vertex to tetrahedron
    e=t4_e2(t);
    v2e = t4util_v_to_t4e_map (size(v,1),e);% Map vertex to edge
    
    % Figure out the  vertex layer numbers which are the guides to coarsening
    if (isempty(vlayer))
        % % Extract the boundary
        f = t4util_boundary(t);
        vlayer =zeros(size(v,1),1,'int16');
        if (surface_coarsening)
            i=setdiff((1:size(v,1))',reshape(f, [],1));
            vlayer(i) = 1;% Mark all vertices in the interior (layer 1) 
            vlayer(f) = 2;% Mark all vertices on the boundary (layer 2) 
            if (preserve_thin) % thin features should be preserved, 
            % mark them as interior (layer 1) so the algorithm does not touch them
                oldvlayer=vlayer;
                for j=reshape(f, 1,[])
                    % is  any vertex connected to j  in the interior?
                    connected =false;
                    for k=1:length(v2e{j})
                        if (oldvlayer(e(v2e{j}(k),1)) == 1) 
                            connected = true; break;
                        end
                        if (oldvlayer(e(v2e{j}(k),2)) == 1) 
                            connected = true; break;
                        end
                    end
                    if (~connected )
                        vlayer(j)=1; % mark it as interior to prevent coarsening
                    end
                end
            end
        else
            vlayer(f) = 1;% Mark all vertices on the boundary (layer 1) 
            % The remaining (interior) vertices will be numbered below 
            % by  distance from the boundary 
        end
    end
    if isempty(find(vlayer~=0))
        return; % I have no information about the layers: bail out
    end
    
    % Compute the vertex layer numbers for vertices for which they have not 
    % been supplied (which is marked by 0).
    currvlayer =1;
    while true
        %mark layer next to current boundary
        oldvlayer=vlayer;
        for i=1:size(vlayer,1)
            if oldvlayer(i)==currvlayer
                for j=1:length(v2t{i})
                    k=v2t{i}(j);
                    for m=1:4
                        if (oldvlayer(t(k,m))==0)
                            vlayer(t(k,m))=currvlayer+1;
                        end
                    end
                end
            end
        end
        if (length(find(vlayer==0))==0) || (norm(double(vlayer-oldvlayer))==0)
            break;
        end
        currvlayer =currvlayer +1;
    end
    currvlayer =currvlayer +1;
    %     Compute the desired mesh size at a given layer
    Layer_desired_es=zeros(currvlayer,1);
    Layer_desired_es(1)=desired_ts;
    for layer =2:currvlayer
        s=0;
        for r=1:currvlayer
            s=s+stretch^r;
            if (s>= layer-1), break,end
        end
        r=r-1;
        Layer_desired_es(layer)=desired_ts*stretch^r;
    end
    %     Layer_desired_es(3:end+1) =(Layer_desired_es(1:end-1)+Layer_desired_es(2:end))/2;
    
    % Initialize edge lengths, edge vertex counts
    es =zeros(size(e,1),1);
    elayer =zeros(size(e,1),1);
    es =edge_length ((1:size(e,1))');
    elayer =min([vlayer(e(:,1)),vlayer(e(:,2))],[],2);
    minne = 200;    maxne = 400;% 36
    
    
    
    Faile =  [];
    pass =1;
    while true
        [availe,currvlayer]=availelist (elayer,currvlayer,minne,maxne,es,nblayer,Layer_desired_es,Faile);
        %disp(['Coarsening: pass ' num2str(pass) ', currvlayer =' num2str(currvlayer)])
        fprintf(1,[num2str(currvlayer) '.']);
        if (mod(pass,20)==0),fprintf(1,['\n']);end
        if (length(availe)==0), break; end % Done. Hallelujah!
        %         pb=progressbar(['Coarsening: pass ' num2str(pass) ', ' num2str(length(availe)) ' edges, currvlayer =' num2str(currvlayer)]);
        Work = 0;
        while true
            eL = sortedelist (elayer,es,Layer_desired_es,availe);
            nWork =size(eL,1);
            Change =  ~true; ntries =0;
            for i=1:length(eL)
                ntries =ntries +1;
                [Success] = collapse_e(eL(i));
                if (Success)
                    Change = true;  break; % minne = 100;    maxne = 200;
                end
                Faile(end+1) =eL(i);
            end
            if (ntries==length(eL)), end% minne =2*minne; maxne =2*maxne;
            if  (~Change)
                %                 try,   pb=progressbar(1.0,pb);delete(pb);  catch,end
                break;
            end
            Work = Work +1; %pb=progressbar(Work/nWork,pb);
        end
        pass = pass +1;
    end
    % % Cleanup
    fprintf(1,['\n'])
    
    [t,v,tmid] =  delete_ent (t,v,tmid);
    
    function [Success] = collapse_e (dei)
        Success = false;
        % if the operation would result in inverted tetrahedra, cancel it
        de =e(dei, [1,2]);
        if (any_negative_volume(t(v2t{de(2)},:),v,de(2),v(de(1),:)))
            de =e(dei, [2,1]);% Try the edge the other way
            if (any_negative_volume(t(v2t{de(2)},:),v,de(2),v(de(1),:)))
                return;% return failure indicator Success
            end
        end
        vi1 =de(1); vi2 =de(2);
        mtl =unique([v2t{vi1},v2t{vi2}]);
        % switch the references to the replaced vertex
        lt=t(mtl,:);
        for k=1:length(mtl)
            if (lt(k,1)==vi2); lt(k,1)=vi1; end;
            if (lt(k,2)==vi2); lt(k,2)=vi1; end;
            if (lt(k,3)==vi2); lt(k,3)=vi1; end;
            if (lt(k,4)==vi2); lt(k,4)=vi1; end;
        end
        t(mtl,:) =lt;
        mel =unique([v2e{vi1},v2e{vi2}]);
        % switch the references to the replaced vertex
        le=e(mel,:);
        for k=1:length(mel)
            if (le(k,1)==vi2), le(k,1) =vi1; end; if (le(k,2)==vi2), le(k,2) =vi1; end;
        end
        e(mel,:)=le;
        dtl =intersect(v2t{vi1},v2t{vi2});
        vl = unique (t(dtl,:));% vertices incident on the collapsed tetrahedra
        for i=1:length(vl) % Delete the collapsed tetrahedra
            tf = ~(ismember(v2t{vl(i)},dtl));           v2t{vl(i)} = v2t{vl(i)}(tf);
        end
        % delete edges which are merged by the collapse
        del =v2e{vi2};
        for k =1:length(del)
            i=del(k);
            if (e(i,1)==vi2),ov=e(i,2);
            else ov =e(i,1); end
            if (~isempty(find(vl==ov)))
                e(i,:) =0;% mark as deleted
                es(i) =inf;% indicate the deleted edge
                elayer(i) =0;% indicate the deleted edge
            end
        end
        t(dtl,:) = 0;% Mark deleted tetrahedra
        e(dei,:) = 0;% Mark deleted edge
        tf = ~(ismember(mtl,dtl));
        v2t{vi1}= mtl(tf);
        v2t{vi2}= [];% this vertex is gone
        tf = ~(ismember(mel,dei));
        v2e{vi1}= mel(tf);
        v2e{vi2}= [];% this vertex is gone
        %     v(vi1,:) = nv1;% new vertex location
        v(vi2,:) = [Inf,Inf,Inf];% Indicate invalid vertex
        % update edge lengths
        for k=1:length(v2e{vi1})
            i=v2e{vi1}(k);
            if (e(i,1)==0),
                es(i) =inf;% indicate the deleted edge
                elayer(i) =0;
            else
                es(i) =edge_length(i);
                elayer(i) =min(vlayer(e(i,:)));
            end
        end
        es(dei) =inf;% indicate the deleted edge
        elayer(dei) =0; % indicate the deleted edge
        Success = true;
    end
    
    
    function eLengths = edge_length (edgeNumbers)
        p=(v(e(edgeNumbers,1),:) - v(e(edgeNumbers,2),:));
        eLengths=sqrt(sum(p.^2,2));
        if (~isempty(vertex_weight))
            vw=max([vertex_weight(e(edgeNumbers,1)),vertex_weight(e(edgeNumbers,2))],[],2);
            eLengths =vw.*eLengths;
        end
    end

end

function eList = sortedelist (elayer,es,Layer_desired_es,availe)
    eList =0*availe;
    n=0;
    for i11=1:length(availe)
        k11=availe(i11);
        if (elayer(k11) >0)
            if  (es(k11)<Layer_desired_es(elayer(k11)))
                n=n+1; eList(n) =k11;
            end
        end
    end
    eList =eList(1:n);
end

function [availe,currvlayer]=availelist (elayer,currvlayer,minnt,maxnt,es,nblayer,Layer_desired_es,Faile)
    eList= [];
    for layer =currvlayer:-1:nblayer+1% This can be changed to allow for more or less coarsening
        availe =setdiff(find(elayer>=layer),Faile);
        eList = sortedelist (elayer,es,Layer_desired_es,availe);
        if (length(eList) >=minnt), break;end
    end
    availe =eList;
    availe =availe (1:min( [length(availe),maxnt]));
    currvlayer =layer;
end


function e=t4_e2(t)
    ec = [1, 2; 2, 3; 3, 1; 4, 1; 4, 2; 4, 3];
    e= [t(:,ec(1,:));t(:,ec(2,:));t(:,ec(3,:));t(:,ec(4,:));t(:,ec(5,:));t(:,ec(6,:));];
    e=sort(e,2);
    [ ignore,ix] =sort(e(:,1));
    ue=e(ix,:);
    e= ue;
    i = 1;
    n=1;
    while n<=size(e,1)
        c =ue(n,1);
        m=n+1;
        while m<=size(e,1)
            if (ue(m,1)~=c), break; end
            m=m+1;
        end
        %  s=unique(ue(n:m-1,2));% the transcription below is a lot faster
        TMPa = ue(n:m-1,2);
        TMPb = sort(TMPa);
        TMPdb = diff(TMPb);
        TMPd = TMPdb ~= 0;
        TMPd(numel(TMPa),1) = true;
        s = TMPb(TMPd);
        % the above lines transcribe unique
        ls =length(s);
        e(i:i+ls-1,1) =c;
        e(i:i+ls-1,2) =s;
        i=i+ls;
        n=m;
    end
    e=e(1:i-1,:);
end



function vol = tetv(X)
    D=X(1,:);
    A=X(2,:)-D; B =X(3,:)-D; C =X(4,:)-D;
    vol =(-A(3)*B(2)+A(2)*B(3))*C(1)+...
        (A(3)*B(1)-A(1)*B(3))*C(2) +...
        (-A(2)*B(1)+A(1)*B(2))*C(3);
end

function vol = tetv1(X)
    vol =det(diff(X));
    %     vol1 =tetv(X);
    %         if abs(vol - vol1)>10*eps,
    %             vol, tetv(X),
    %             disp('Different numbers')
    %         end
    %         if (vol * vol1)<0,
    %             vol, tetv(X),
    %             disp('Opposite signs')
    %         end
end

function bul =any_negative_volume(t,v,vi2,v1)
    for iS1 =1:size(t,1)
        lt=t(iS1,:);
        lv =v(lt,:);
        if (lt(1)==vi2); lv(1,:)=v1; end;
        if (lt(2)==vi2); lv(2,:)=v1; end;
        if (lt(3)==vi2); lv(3,:)=v1; end;
        if (lt(4)==vi2); lv(4,:)=v1; end;
        if (tetv1(lv)<0),bul = true; return; end
    end
    bul= false;
end

function [t,v,tmid] =  delete_ent (t,v,tmid)
% Volumes = 0*double(t(:,1));
%     for iS1 =1:size(t,1)
%         c=t(iS1,:);
%         if (~any(c<1))
%         Volumes (iS1) =tetvol(v(c,:));
%         end
%     end
% if (any(Volumes<0))
% disp( 'Negative volume' )
% end
    
    nn =zeros(size(v,1),1);
    nv =0*v;
    k=1;
    for i=1:size(v,1)
        if (v(i,1)~=inf)
            nn(i) =k;
            nv(k,:) =v(i,:);
            k=k+1;
        end
    end
    nnv=k-1;
    v=nv(1:nnv,:);
    % delete connectivities of collapsed tetrahedra, and renumber nodes
    nt=0*t;
    ntmid=0*tmid;
    k=1;
    for i=1:size(t,1)
        if (t(i,:)~=0)% not deleted
            if (~isempty(find(nn(t(i,:))==0)))
                %                 error('Referring to deleted vertex')
                t(i,:)=0;
            else
                j =nn(t(i,:));
                try
                    nt(k,:)=j;
                    ntmid(k)=tmid(i);
                    k=k+1;
                catch
                end
            end
        end
    end
    t=nt(1:k-1,:);
    % delete unconnected vertices
    uv = unique([t(:,1);t(:,2);t(:,3);t(:,4)]);
    if (length(uv)~=size(v,1))
        nn =zeros(size(v,1),1);
        nn(uv) =(1:length(uv))';
        nv =0*v;
        for i=1:size(v,1)
            if (nn(i)~=0)
                nv(nn(i),:) =v(i,:);
            end
        end
        v=nv(1:length(uv),:);
        for i=1:size(t,1)
            t(i,:) =nn(t(i,:));
        end
    end
    tmid=ntmid(1:k-1);
end

% function vol = tetvol(X)
%         vol = det([X(2,:)-X(1,:);X(3,:)-X(1,:);X(4,:)-X(1,:)])/6;
%     end
