function d= extract_displacement_from_abaqus_dat(File,Tagline,nnodes)
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

