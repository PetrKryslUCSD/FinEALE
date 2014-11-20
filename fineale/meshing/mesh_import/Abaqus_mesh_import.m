function [fens,fes]=Abaqus_mesh_import(filename,xyzscale,saveReadData)
% Import mesh from the Abaqus mesh file.
%
% function [fens,fes]=Abaqus_mesh_import(filename,xyzscale,saveReadData) 
%
% Import tetrahedral (4- and 10-node) ABAQUS Mesh (.INP).
% Limitations: only the *NODE and *ELEMENT  sections are read. Only 4-node and 10-node
% tetrahedra  are handled.
%
% Arguments:
% filename= filename,
% xyzscale= scale the coordinates by this number when they are read in,
% saveReadData= should the data of the saved in a format that allows 
%           for speedy reading of the same data later on?
% 
% Output:
% fens= finite element node set
% fes = finite element set
%
%
% See also: mesh_import
%    
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
        
        chunk=1000;
        
        nnode=0;
        node=zeros(chunk,4);
        temp = '';
        while ~strncmpi(temp,'*NODE',5)
            temp = fgetl(fid);   if ~ischar(temp), break, end
        end
        if (~strncmpi(temp,'*NODE',5))
            error
        end
        
        temp = '';
        while true
            nnode=nnode+1;
            temp = fgetl(fid);
            if (strncmpi(temp,'*',1))
                nnode=nnode-1;
                break
            end
            A = sscanf(temp, '%g,%g,%g,%g');
            node(nnode,:)=A;
            if (size(node,1)<=nnode)
                node(size(node,1)+1:size(node,1)+chunk,:)=zeros(chunk,4);
            end
        end
        
        % Now find and process all *ELEMENT blocks
        More_data=true;
        elsets = {};
        while More_data
            % Find the next block
            while ~strncmpi(temp,'*ELEMENT',8)
                temp = fgetl(fid);
                if ~ischar(temp),
                    More_data=~true;
                    break,
                end
            end
            [Valid, Type, Elset]=Parse_element_line(temp);
            if (Valid)
                % Valid element type
                elsets{end+1}.Elset=Elset;
                elsets{end}.nelem=0;
                elsets{end}.elem= [];
                elsets{end}.ennod= Type;
                temp = '';
                while true
                    elsets{end}.nelem=elsets{end}.nelem+1;
                    temp = fgetl(fid);
                    if (isempty(elsets{end}.elem))
                        All =sscanf(temp, '%d,',inf);
                        if (elsets{end}.ennod ~=length(All)-1)
                            error(' Wrong number of data items');
                        end
                        elsets{end}.elem=zeros(chunk,elsets{end}.ennod+1);
                    end
                    if (strncmpi(temp,'*',1))
                        elsets{end}.nelem=elsets{end}.nelem-1;
                        break
                    end
                    A = sscanf(temp, '%d,',inf);
                    elsets{end}.elem(elsets{end}.nelem,:)=A';
                    if (size(elsets{end}.elem,1)<=elsets{end}.nelem)
                        elsets{end}.elem(size(elsets{end}.elem,1)+1:size(elsets{end}.elem,1)+chunk,:)=zeros(chunk,elsets{end}.ennod+1);
                    end
                end
            end
        end
        
        fclose(fid);
        
        if (norm((1:nnode)'-node(1:nnode,1)))
            error
        end
        
        % Process output arguments
        % Extract coordinates
        xyz=node(node(1:nnode,1),2:4);
        xyz =xyzscale*xyz;
        % Cleanup element connectivities
        for j=1:length(elsets)
            elsets{j}.conn=elsets{j}.elem(1:elsets{j}.nelem,2:end);
        end
    end
    
    % Create output arguments. First the nodes
    clear fens
    fens=fenode_set(struct('xyz',xyz));
    
    % Now the geometric cells for each element set
    clear fes
    for j=1:length(elsets)
        switch elsets{j}.ennod
            case 4
                Constructor =@fe_set_T4;
            case 10
                Constructor =@fe_set_T10;
            otherwise
                error('Unable to open specified file.');
        end
        fes{j}= Constructor(struct('conn',elsets{j}.conn(:,1:elsets{j}.ennod),'label',j));
    end
    % If there was just one element set, return a single fe set object
    if (length(elsets)==1)
        fes=fes{1};
    end
    
    % Should the data that we read be saved?
    if (saveReadData && (Doread))
        save(saveFile,'xyz','elsets')
    end
end
function [Valid, Type, Elset]=Parse_element_line(String)
    Valid=false; Type= []; Elset=[];
    if ~ischar(String), return; end
    tokens= all_tokens(String, ',');
    if (strcmpi(tokens{1},'*ELEMENT'))
        if (length(tokens)>=2)
            tok1= all_tokens(tokens{2}, '=');
            if (strcmpi(tok1{1},'TYPE'))
                if (strcmpi(tok1{2},'C3D4'))
                    Type=4;
                elseif (strcmpi(tok1{2},'C3D10'))
                    Type=10;
                end
            end
        end
        if (length(tokens)>=3)
            tok1= all_tokens(tokens{3}, '=');
            if (strcmpi(tok1{1},'ELSET'))
                Elset=tok1{2};
            end
        end
        Valid=true;
    end
end
function tokens= all_tokens(String, Separator)
    tokens={};
    while ~isempty( String )
        [tokens{end+1}, String]=strtok(regexprep(String,'\s*',''),Separator);
    end
end
