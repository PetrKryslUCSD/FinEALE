function mlar(theFile,varargin)
% Matlab plain-text archiver.
%
% function mlar(theFile,varargin)
% 
% Invoking with two or more arguments will treat the second, third, and so
% on argument as file names to be archived. The result is stored in the
% file named  as given in the argument theFile.
% 
% Invoking with only one argument runs the archiver on the supplied file to
% extract the stored files from it.
%
% Run without arguments extracts all archived files from the present file.
% 
% Example:
% mlar t.m j23.jpg
% Stores the file j23.jpg into the archive t.m
%
% mlar Test.m
% Extracts all the files stored in Test.m
%
% mlar
% Extracts all the files stored in mlar.m. 
% 
% Note that the name of this file may be different from the name of the
% function. This would be the case when this file was written by
% make_patch(). Files would have been attached to the end of mlar.m, which
% would then have been named something else (for instance p1.m). The stored
% files then could be extracted from p1.m by running it without arguments.
%
%
% See also: make_patch
%
% Copyright(C) 2007, 2008, 2011 Petr Krysl
    if (nargin<2)%
        if (nargin==0)%
            theFile =[mfilename,'.m'];
        end
        fid=fopen(theFile,'r');
        if (fid==-1)
            warning (['Could not open ' theFile])
            return
        end
        mlar_Read(fid);
        fid=fclose(fid);
    else
        for i= 2:nargin
            fid=fopen(theFile,'a');
            if (fid==-1)
                warning (['Could not open ' theFile])
            else
                mlar_Write(fid,varargin{i-1});
                fid=fclose(fid);
            end
        end
    end
end

function mlar_Write(fid,theFile)
    [isbin,type] = isbinary(theFile);

    if ~isempty(isbin)
        ifid=fopen(theFile);
        if (ifid==-1)
            warning (['Could not open ' theFile])
            return
        else
            disp (['Archiving ' theFile])
        end
        c=fread(ifid);
        ifid=fclose(ifid);
        theRealSize=prod(size((c)));
        if ~isbin
            c=char(c');
            c=strrep(c,char(13),'');
            if c(end)==10, c=c(1:end-1); end
            c=[char(37),strrep(c,char(10),[10 char(37)]),char(10)];
            theType = 'text';
        else
            c=sprintf('%02x\n',c);
            c(3:3:end)='';
            nChar=70;
            m=floor(length(c)/nChar);
            tmp='';
            tmp=reshape(c(1:m*nChar),nChar,m);
            tmp=tmp';
            if mod(length(c),nChar)
                tmp=strvcat(tmp,c(m*nChar+1:end));
            end
            tmp(:,2:end+1)=tmp;
            tmp(:,1)=char(37);% insert a %
            tmp(:,end+1) = 10;
            tmp(1,end)=10;
            c=tmp;
            theType = 'binary';
        end
        theSize=prod(size((c)));
        if ~isbin
            fprintf(fid,'\n');
            fprintf(fid,'%%MLAR ASCII begin:%s\n',theFile);
            fwrite(fid,c);
            fprintf(fid,'%%MLAR ASCII end:%s\n',theFile);
        else
            fprintf(fid,'\n');
            fprintf(fid,'%%MLAR BIN begin:%s\n',theFile);
            fprintf(fid,'%s',c');
            fprintf(fid,'%%MLAR BIN end:%s\n',theFile);
        end

    end
end

function [isbin,type]=isbinary(file)
    isbin=0;
    type='';
    fid=fopen(file);
    if fid~=-1
        c=fread(fid,1000);
        if any((c>=0) & (c<9)), isbin = 1; end
        if any((c>=14) & (c<32)), isbin = 1; end
        fclose(fid);
        if isbin
            type='binary';
        else
            type='text';
        end
    else
        isbin=[];
        type='';
    end
end

function  mlar_Read(fid)
    Reading_ASCII = ~true;
    Reading_BIN = ~true;
    while true
        tline = fgetl(fid);
        if ~ischar(tline),   break,   end
        if (Reading_ASCII)
            if (length(tline)>=1) && (tline(1) == char(37))
                if (length(tline)>=16) && (strcmp(tline(2:16),'MLAR ASCII end:'))
                    Reading_ASCII = ~true;
                    ofid=fclose(ofid);
                else
                    fprintf(ofid, '%s\n', tline(2:end));
                end
            else
                warning(['Line does not start with a %: ''' tline '''']);
                fprintf(ofid, '%s\n', ['****** ' tline(2:end)]);
            end
        elseif (Reading_BIN)
            if (length(tline)>=1) && (tline(1) == char(37))
                if (length(tline)>=14) && (strcmp(tline(2:14),'MLAR BIN end:'))
                    Reading_BIN = ~true;
                    ofid=fclose(ofid);
                else
                    fwrite(ofid,h2b(tline(2:end)));
                end
            else
                warning(['Line does not start with a %: ''' tline '''']);
                fprintf(ofid, '%s\n', ['****** ' tline(2:end)]);
            end
        else
            if  (length(tline)>=18) && (strcmp(tline(2:18),'MLAR ASCII begin:'))
                [pathstr, name, ext, versn] = fileparts(tline(19:end)) ;
                if (isempty(pathstr) ||  (pathstr(1)~=filesep))
                    pathstr = ['.',filesep,pathstr];
                end
                if (exist(pathstr,'dir')~=7)
                    [s, mess, messid] = mkdir(pathstr);
                end
                ofid=fopen(tline(19:end),'w');
                if (ofid==-1)
                    warning (['Could not open ' tline(19:end)])
                    Reading_ASCII =  ~true;
                else
                    disp (['Extracting ' tline(19:end)])
                    Reading_ASCII = true;
                end
            elseif (length(tline)>=16) && (strcmp(tline(2:16),'MLAR BIN begin:'))
                [pathstr, name, ext, versn] = fileparts(tline(17:end)) ;
                if (pathstr(1)~=filesep)
                    pathstr = ['.',filesep,pathstr];
                end
                if (exist(pathstr,'dir')~=7)
                    [s, mess, messid] = mkdir(pathstr);
                end
                ofid=fopen(tline(17:end),'w');
                if (ofid==-1)
                    warning (['Could not open ' tline(17:end)])
                    Reading_BIN =  ~true;
                else
                    disp (['Extracting ' tline(17:end)])
                    Reading_BIN = true;
                end
            end
        end
    end
end

function b=h2b(h)
    c = zeros(1, 256);
    c(abs('0'):abs('9')) = 0:9;
    c(abs('a'):abs('f')) = 10:15;
    hex = double(h);
    b = 16*c(hex(1:2:end)) + c(hex(2:2:end));
end


