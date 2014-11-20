function sy= remove_bias(t,y)
    % Try to remove  linear drift (bias) from a  signal.
    %
    %  function sy= remove_bias(t,y)

     [m,n]=size(t );
    t=t/mean(t);
    t=reshape(t,length(t),1);
    y=reshape(y,length(y),1);
    %     if (length(ts)>2)
    %         B=[ts,0*ts+1];
    %         p = (B'*B)\(B'*ys);
    %     else 
        B=[t,0*t+1];
        p = (B'*B)\(B'*y);
        %     end 
    sy=y-p(1)*t.^1-p(2);
    %     sy= remove_bias_from_peaks(t,sy); 
end

function sy= remove_bias_from_peaks(t,y)
    [m,n]=size(t );
    t=t/mean(t);
    t=reshape(t,length(t),1);
    y=reshape(y,length(y),1);
    [ppeaks,ppeakts,npeaks,npeakts]=  find_peaks(t,y);
    ts=[ppeakts,npeakts];
    ys=[ppeaks,npeaks];
    ts=reshape(ts,length(ts),1);
    ys=reshape(ys,length(ys),1);
    if (length(ts)>4)
        B=[ts,0*ts+1];
        p = (B'*B)\(B'*ys);
    else
        B=[t,0*t+1];
        p = (B'*B)\(B'*y);
    end
    sy=y-p(1)*t.^1-p(2);
end


 %         B=[0*t+1];
   %         B=[t.^3,t.^2,t,0*t+1];
    %         B=[t.^2,t,0*t+1];
    %         sy=y-p(1)*t.^3-p(2)*t.^2-p(3)*t-p(4);
    % sy=y-p(1)*t.^2-p(2)*t-p(3);
       %         sy=y-p(1);
 
function [ppeaks,ppeakts,npeaks,npeakts]= find_peaks(t,y)
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