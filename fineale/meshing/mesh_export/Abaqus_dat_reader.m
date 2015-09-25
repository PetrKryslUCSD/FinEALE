classdef  Abaqus_dat_reader  < handle
    %     Simple class to read an abaqus DAT file.
    %
   
    
    methods % constructor
        
        function self = Abaqus_dat_reader(Parameters)
            if nargin <1
                Parameters = struct([]);
            end
            
        end
        
    end
    
    
    methods % concrete methods
        
        function A= extract_energy_from_abaqus_dat(self,File,Tagline)
            A=[];
            try
                fid=fopen(File,'r');
                while true
                    temp = fgetl(fid);   if ~ischar(temp), break, end
                    temp=deblank(fliplr(deblank(temp(end:-1:1))));
                    if (strcmpi(temp,Tagline))
                        for j=1:3
                            temp = fgetl(fid);
                        end
                        temp = fgetl(fid);
                        A = sscanf(temp(40:end), '%g');
                        fclose(fid);
                        return;
                    end
                end
            catch,fclose(fid);end
        end
        
        function d= extract_displacement_from_abaqus_dat(self,File,Tagline,nnodes)
            if (~exist('nnodes','var') ),nnodes=1;end
            d=zeros(nnodes,3);
            try
                fid=fopen(File,'r');
                while true
                    temp = fgetl(fid);   if ~ischar(temp), break, end
                    temp=deblank(fliplr(deblank(temp(end:-1:1))));
                    if (strcmpi(temp,Tagline))
                        for j=1:4
                            temp = fgetl(fid);
                        end
                        for j=1:nnodes
                            temp = fgetl(fid);
                            A = sscanf(temp, '%g %g %g %g');
                            d(j,:)=A(2:end);
                        end
                        fclose(fid);
                        return;
                    end
                end
            catch,fclose(fid);end
        end
        
        function d= extract_buckling_from_abaqus_dat(self,File,Tagline,neigv)
            if (~exist('neigv','var') ),neigv=1;end
            d=zeros(neigv,1);
            try
                fid=fopen(File,'r');
                while true
                    temp = fgetl(fid);   if ~ischar(temp), break, end
                    temp=deblank(fliplr(deblank(temp(end:-1:1))));
                    if (strcmpi(temp,Tagline))
                        for j=1:2
                            temp = fgetl(fid);
                        end
                        for j=1:neigv
                            temp = fgetl(fid);
                            A = sscanf(temp, '%g %g');
                            d(j,:)=A(2:end);
                        end
                        fclose(fid);
                        return;
                    end
                end
            catch,fclose(fid);end
        end
        
        
        
    end
    
end


