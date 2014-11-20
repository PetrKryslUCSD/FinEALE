function draw_grid_lines(self, options)
% Draw the grid lines of the global Cartesian coordinate system.
% THIS ISN'T WORKING YET!
%
% function draw_grid_lines(self, options)
%
% Draw the grid lines of the global Cartesian coordinate system.
% 
% options= optional structure, with fields
% length = length of the axes' arrows
% 
    plot3(get(gca,'XLim'),[0 0],[0 0],'k');
    plot3([0 0],[0 0],get(gca,'ZLim'),'k');
    plot3([0 0],get(gca,'YLim'),[0 0],'k');

    % GET TICKS
    X=get(gca,'Xtick');
    Y=get(gca,'Ytick');
    Z=get(gca,'Ztick');

    % GET LABELS
    XL=get(gca,'XtickLabel');
    YL=get(gca,'YtickLabel');
    ZL=get(gca,'ZtickLabel');

    % REMOVE TICKS
    set(gca,'Xtick',[]);
    set(gca,'Ytick',[]);
    set(gca,'Ztick',[]);

    % GET OFFSETS
    Xoff=diff(get(gca,'XLim'))./30;
    Yoff=diff(get(gca,'YLim'))./30;
    Zoff=diff(get(gca,'ZLim'))./30;

    % DRAW TICKS
    %%%%%%% THIS COULD BE VECTORiZeD %%%%%%%
    for i=1:length(X)
        plot3([X(i) X(i)],[0 0],[-Zoff Zoff],'k');
    end;
    for i=1:length(Y)
        plot3([-Xoff Xoff],[Y(i) Y(i)],[0 0],'k');
    end;
    for i=1:length(Z)
        plot3([-Xoff Xoff],[0 0],[Z(i) Z(i)],'k');
    end;

    % DRAW LABELS
    text(X,zeros(size(X)),zeros(size(X))-3.*Zoff,XL);
    text(zeros(size(Y))-3.*Xoff,Y,zeros(size(Y)),YL);
    text(zeros(size(Z))-3.*Xoff,zeros(size(Z)),Z,ZL);
end

