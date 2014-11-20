function n=frame_name(p,f,i,s)
    if (strcmp(class(i),'char'))
        n=[p filesep f i '.' s];
    else
        if i <10
            n=[p filesep f '-0000' num2str(i) '.' s];
        elseif i<100
            n=[p filesep f '-000' num2str(i) '.' s];
        elseif i<1000
            n=[p filesep f '-00' num2str(i) '.' s];
        elseif i<10000
            n=[p filesep f '-0' num2str(i) '.' s];
        else
            n=[p filesep f '-' num2str(i) '.' s];
        end
    end
end
