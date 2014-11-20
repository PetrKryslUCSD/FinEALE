function [fens,fes]=Ansys_mesh_import(filename,xyzscale,saveReadData)
    % Import mesh from the Ansys mesh file.
    %
    % function [fens,fes]=Ansys_mesh_import(filename,xyzscale,saveReadData)
    %
    % Import tetrahedral (4- and 10-node) Ansys Mesh (.cdb).
    %
    % Limitations:
    % 1. Only the N and EN commands are processed.
    % 2. Only 4-node and 10-node tetrahedra  are handled.
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
    % Example:
    % [fens,fes]=Ansys_mesh_import('cfin.cdb');%Note: we know there are three element sets in this file
    % gv=drawmesh({fens,fes{1}},'fes','facecolor','r')
    % gv=drawmesh({fens,fes{2}},'fes','facecolor','g','gv',gv)
    % gv=drawmesh({fens,fes{3}},'fes','facecolor','b','gv',gv)
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
        elsets = {};
        
        while true
            temp = fgetl(fid);
            if ~ischar(temp),
                break,
            end
            if (strncmpi(temp,'N,R',3))%  The N command
                revision=sscanf(temp, 'N,R%g');% Revision number
                if (revision == 5.0)
                    % N,R5.x,Type,Node,Solid,Parm,Val1,Val2,Val3
                    % Type =LOC...Read in coordinates
                    %     ANG...Read in rotation angles
                    % Node= node number
                    %  Solid= solid model reference key, not present for ANG
                    %  Parm= line parameter value, not present for ANG
                    % Val1,Val2,Val3P = values to be read
                    tokens= all_tokens(temp, ',');
                    if (strcmpi(tokens{3},'LOC'))
                        ID =str2num(tokens{4}); XYZ=[str2num(tokens{5}),str2num(tokens{6}),str2num(tokens{7})];;
                    else
                        error(['Could not parse ' temp]);
                    end
                else
                    error(['Revision not recognized:' temp]);
                end
                %
                nnode=nnode+1;
                node(nnode,:)=[ID, XYZ];
                if (size(node,1)<=nnode)
                    node(size(node,1)+1:size(node,1)+chunk,:)=zeros(chunk,4);
                end
            elseif (strncmpi(temp,'EN,R',4))%  The EN command
                revision=sscanf(temp, 'EN,R%g');% Revision number
                if (revision == 5.0)||(revision == 5.1)
                    % EN,R5.x,Type,Numn,Val1...Val8
                    % Type =ATTR...Read in attributes
                    %     NODE...Read in connectivity
                    % For  Type =ATTR : NUMN,MAT,TYPE,REAL,NUMELEMENT,ESYS, Solid_Ref,Birth/Death,Pexclude
                    tokens= all_tokens(temp, ',');
                    if (strcmpi(tokens{3},'ATTR'))
                        NUMN=str2num(tokens{4});
                        MAT=str2num(tokens{5});
                        TYPE=str2num(tokens{6});
                        NUMELEMENT=str2num(tokens{8});
                        elsetn=MAT; %  We are assuming here something about assignment of elements to element sets
                        %                         If the element set number is the type, change the above line.
                        resize_elsets(elsetn,NUMELEMENT);
                        %             elsets{elsetn}.conn(NUMELEMENT,:)= [];
                        if (elsets{elsetn}.NUMN ~= NUMN)...
                                || (elsets{elsetn}.MAT ~= MAT)...
                                || (elsets{elsetn}.TYPE ~= TYPE)
                            error(['Element set mismatch ' temp]);
                        end
                        % This command should be followed by the EN command with Type=  NODE
                        n=0;;
                        conn=zeros(1,NUMN);
                        while n<NUMN
                            temp = fgetl(fid);
                            if ~ischar(temp),
                                break,
                            end
                            if (strncmpi(temp,'EN,R',4))%  The EN command
                                tokens= all_tokens(temp, ',');
                                if (strcmpi(tokens{3},'NODE'))
                                    for j=1:min([NUMN-n,8])
                                        conn(1,n+j)=str2num(tokens{3+j});
                                    end
                                    n=n+min([NUMN-n,8]);
                                else
                                    error(['Expected EN command with Type=  NODE, got ' temp]);
                                end
                            else
                                error(['Expected EN command with Type=  NODE, got ' temp]);
                            end
                        end
                        elsets{elsetn}.conn(NUMELEMENT,:)=conn;
                    else
                        error(['Expected EN command with Type=  ATTR, got ' temp]);
                    end
                else
                    error(['Revision not recognized:' temp]);
                end
            end
        end
        
        fclose(fid);
        
        %Rearrange the nodes so that they are in the order of their numbers
        %Also, trim out the rows not assigned
        node=node(node(1:nnode,1),:);
        if (norm((1:nnode)'-node(1:nnode,1)))
             error(['Node numbering inconsistent' ]);
        end
        
        % Process output arguments
        % Extract coordinates
        xyz=node(node(1:nnode,1),2:4);
        xyz =xyzscale*xyz;
        % Cleanup element connectivities
        for j=1:length(elsets)
            ix=find(elsets{j}.conn(:,1)~=0);
            elsets{j}.conn=elsets{j}.conn(ix,:);
        end
    end
    
    % Create output arguments. First the nodes
    clear fens
    fens=fenode_set(struct('xyz',xyz));
    
    % Now the geometric cells for each element set
    clear fes
    for j=1:length(elsets)
        switch elsets{j}.NUMN
            case 4
                Constructor =@fe_set_T4;
            case 10
                Constructor =@fe_set_T10;
            otherwise
                error('Unable to open specified file.');
        end
        fes{j}= Constructor(struct('conn',elsets{j}.conn(:,1:elsets{j}.NUMN),'label',j));
    end
    % If there was just one element set, return a single fe set object
    if (length(elsets)==1)
        fes=fes{1};
    end
    
    % Should the data that we read be saved?
    if (saveReadData && (Doread))
        save(saveFile,'xyz','elsets')
    end
    
    function resize_elsets(elsetn,nelem)
        if (isempty(elsets)) || (length(elsets)<elsetn)
            elsets{elsetn}.conn= [];
            elsets{elsetn}.NUMN= NUMN;
            elsets{elsetn}.MAT= MAT;
            elsets{elsetn}.TYPE= TYPE;
        end
        if (size(elsets{elsetn}.conn,1)<=nelem)
            elsets{elsetn}.conn(size(elsets{elsetn}.conn,1)+1:size(elsets{elsetn}.conn,1)+chunk,:)=...
                zeros(chunk,elsets{elsetn}.NUMN);
        end
    end
end


function tokens= all_tokens(String, Separator)
    tokens={};
    while ~isempty( String )
        [tokens{end+1}, String]=strtok(regexprep(String,'\s*',''),Separator);
    end
end
