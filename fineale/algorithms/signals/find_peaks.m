
function [ppeaks,ppeakts,npeaks,npeakts]= find_peaks(t,y)
    % Find peaks in a (harmonic) signal.
    %
    % function [ppeaks,ppeakts,npeaks,npeakts]= find_peaks(t,y)
    %
    %     The function tries to find positive and negative peaks, both the
    %     values and the times at which the peaks occurred.
    ppeaks=[]; npeaks=[];
    ppeakts=[]; npeakts=[];
    dy=diff(y);
    for j3=1:length(dy)-1
        if (dy(j3)*dy(j3+1)<0)
            if (dy(j3)<0)
                npeaks(end+1)=y(j3+1);
                npeakts(end+1)=t(j3+1);
            else
                ppeaks(end+1)=y(j3+1);
                ppeakts(end+1)=t(j3+1);
            end
        end
    end
end