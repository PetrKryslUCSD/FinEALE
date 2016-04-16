% Import tetrahedral (4- and 10-node) NASTRAN mesh.
%
% function [fens,gcells]=nastran_mesh_import(filename,xyzscale,saveReadData)
%
% Limitations:
% 1. only the GRID and CTETRA  sections are read. 
% 2. Only 4-node and 10-node tetrahedra  are handled.
% 3.  The file needs to be free-form (data separated by commas).
function [fens,gcells]=nastran_mesh_import(filename,xyzscale,saveReadData)
    
    if (~exist('xyzscale','var'))
        xyzscale = 1.0;
    end
    
    if (~exist('saveReadData','var'))
        saveReadData = ~true;
    end
    
    [pathstr, name, ext] = fileparts(filename) ;
    if (isempty(pathstr)),pathstr ='.';end
    saveFile = [pathstr filesep name ,'_xyzconn' '.mat'];
    
    Doread = true;
    if (saveReadData)
        if (exist(saveFile,'file'))
            try,
                load(saveFile);
                ennod =size(conn,2);
                Doread= false;
            catch,end;
        end
    end
    if (Doread)
        fid=fopen(filename,'r');
        
        if fid == -1,
            error('Unable to open specified file.');
        end
        
        chunk=5000;
        
        nnode=0;
        node=zeros(chunk,4);
        nelem=0;
        elem=zeros(chunk,13);
        ennod=[];
        temp = '';
        while true
            temp = fgetl(fid);
            if ~ischar(temp),
                More_data=~true;
                break,
            end
            temp=strtrim(temp);
            if (strncmpi(upper(temp),'GRID',4))
                % Template:
                %                 GRID,1,,-1.32846E-017,3.25378E-033,0.216954
                nnode=nnode+1;
                temp=strrep (temp,',',' ');
                A = sscanf(temp, 'GRID %g %g %g %g');
                node(nnode,:)=A;
                if (size(node,1)<=nnode)
                    node(size(node,1)+1:size(node,1)+chunk,:)=zeros(chunk,4);
                end
            elseif (strncmpi(upper(temp),'CTETRA',6))
                % Template:
                %                 CTETRA,1,3,15,14,12,16,8971,4853,8972,4850,8973,4848
                nelem=nelem+1;
                temp=strrep (temp,',',' ');
                A=Parse_CTETRA_line(temp);
                elem(nelem,1:length(A))=A;
                if (size(elem,1)<=nelem)
                    elem(size(elem,1)+1:size(elem,1)+chunk,:)=zeros(chunk,13);
                end
            elseif (strncmpi(temp,'$',1))
                % Comment line
            elseif (strncmpi(upper(temp),'BEGIN BULK',1))
           end
        end
        
        node=node(1:nnode,:);
        elem=elem(1:nelem,:);
        
        fclose(fid);
        
        if (norm((1:nnode)'-node(1:nnode,1)))
            error('Nodes are not in serial order');
        end
        
        % Process output arguments
        % Extract coordinates
        xyz=node(node(1:nnode,1),2:4);
        xyz =xyzscale*xyz;
        % Cleanup element connectivities
        ennod=unique(elem(:,3));
        if  length(ennod)~=1
            error('This function cannot treat a mixture of element types at this point');
        end
        conn=elem(:,4:4+ennod-1);
        label =elem(:,2);
    end
    
    % Create output arguments. First the nodes
    clear fens
    fens=fenodeset(struct('id',(1:size(xyz,1))','xyz',xyz));
    
    % Now the geometric cells for each element set
    clear gcells
    switch ennod
        case 4
            Constructor =@gcellset_T4;
        case 10
            Constructor =@gcellset_T10;
        otherwise
        error('Unknown element type');
    end
    gcells= Constructor(struct('id',(1:size(conn,1))','conn',conn(:,1:ennod),'label',label));
    
    % Should the data that we read be saved?
    if (saveReadData && (Doread))
        save(saveFile,'xyz','elsets')
    end
end

function A=Parse_CTETRA_line(String)
    A=[];
    if ~ischar(String), return; end
    tokens= all_tokens(String, ' ');
    if (strcmpi(tokens{1},'CTETRA'))
        A(1)=sscanf(tokens{2}, '%g');% Element number
        A(2)=sscanf(tokens{3}, '%g');% Region number
        A(3)=(length(tokens)==7)*4 + (length(tokens)==13)*10;% nodes per element, 4 or 10
        for k=1:A(3)
            A(3+k)=sscanf(tokens{3+k}, '%g');% Node number
        end 
    end
end
function tokens= all_tokens(String, Separator)
    tokens={};
    while ~isempty( String )
        [tokens{end+1}, String]=strtok(String,Separator);
    end
end
