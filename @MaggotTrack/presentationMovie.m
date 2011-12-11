function varargout = presentationMovie(track, varargin)
%@MaggotTrack
%playMovie(track, varargin)
%enter options as pairs, caps matter
%options, with defaults
%
%ptbuffer = 1000;
%delayTime = 0.05;
%inds = 1:length(track.pt);
%startLoc = []; > if startLoc & stopLoc are not empty, we run the movie
%between these two points
%stopLoc = []; >
%pass 'fid', [] to not load images from disk
%'AxesList', ax = list of axes to plot to; specified from top left in cc
% order
%'TitleOptions', {'Color', 'k', 'FontSize', 14};
%'DataPlotOptions', {'Color', 'k', 'LineWidth', 3};
%'DataAxesOptions', {}
%'LocPlotOptions', {'field', 'sloc', 'LineWidth', 2, 'AxesOptions', {}}

ptbuffer = 200;
delayTime = 0.05;

inds = 1:length(track.pt);
track.expt.openDataFile;
fid = track.expt.fid;
AxesList = [];
TitleOptions = {'Color', 'k', 'FontSize', 14};
DataPlotOptions = {'Color', 'k', 'LineWidth', 3};
DataAxesOptions = {};
LocPlotOptions = {'field', 'sloc', 'LineWidth', 2, 'AxesOptions', {}};
aviobj = [];
%vidobj requires matlab 2011 (7.12) or later
vidobj = [];
avirect = [];
varargin = assignApplicable(varargin);
if (isempty(AxesList))
    order = [1 3 4 2];
    for j = 1:4
        AxesList(j) = subplot(2,2,order(j)); %#ok<AGROW>
    end
end
        

pt = [track.pt];
%loc = [pt.loc];
sloc = track.getDerivedQuantity('sloc');
sind = track.getDerivedQuantity('mapptstointerped');
track.calculateDerivedQuantity({'sbodytheta', 'speed', 'vel_dp'});

sstart = sind(1) - ptbuffer;
send = sind(end) + ptbuffer;
if (sstart < 1)
    sstart = 1;
end
if (send > length(sloc))
    send = length(sloc);
end
datafields(track, AxesList([4 2]), sstart:send, TitleOptions, DataPlotOptions, DataAxesOptions);
locPlot(track, sstart:send, AxesList(3), LocPlotOptions{:}, 'TitleOptions', TitleOptions);
handles = [];
hloc = [];



for j = inds
    ts1 = tic();
    axcolor = get(AxesList(1), 'Color');
    pt(j).drawTrackImage(track.expt.camcalinfo,'Axes', AxesList(1), 'fid', fid, varargin{:}); 
    shading(AxesList(1),'interp');
    xl = get(AxesList(1),'XLim');
    yl = get(AxesList(1),'YLim');
    
    hold (AxesList(1),'on');
    sstart = sind(j) - ptbuffer;
    send = sind(j) + ptbuffer;
    if (sstart < 1)
        sstart = 1;
    end
    if (send > length(sloc))
        send = length(sloc);
    end
    
    plot (AxesList(1),sloc(1,sstart:sind(j)), sloc(2,sstart:sind(j)), 'b.-');
    plot (AxesList(1),sloc(1,sind(j)), sloc(2,sind(j)), 'bo', 'MarkerSize', 5);
    set(AxesList(1),'XLim', xl, 'YLim', yl, 'YDir', 'reverse');
    hold (AxesList(1), 'off');
    set(AxesList(1), 'XTick', [], 'YTick', [],'Color',axcolor);
    if (~isempty(track.run))
        t = [];
        
        if (~isempty(track.headSwing) && any([track.headSwing.inds] == sind(j)))
            t = 'headsweep ';
            I = find([track.headSwing.startInd] <= sind(j) & [track.headSwing.endInd] >= sind(j));
            if (~isempty(I))
                if (track.headSwing(I).accepted)
                    t = [t 'accepted ']; %#ok<AGROW>
                else
                    t = [t 'rejected ']; %#ok<AGROW>
                end
            end
        else
            if (track.isrun(sind(j)))
                t = 'run '; 
            end
        end
        %{
        if (~isempty(track.reorientation) && any([track.reorientation.inds] == sind(j)))
            t = [t 'reorientation '];
        end
        %}
        title (AxesList(1),t, TitleOptions{:});
    end
    
    handles = updateCenter(handles, AxesList([4 2]),track, sind(j), sstart, send, DataAxesOptions);
    hloc = locUpdate (track, sind(j), AxesList(3), hloc, LocPlotOptions{:});
    timeleft = delayTime - toc(ts1);
    if (timeleft > 0)
        pause(timeleft);
    else
        pause(0.001);
    end
    if (~isempty(aviobj) || ~isempty(vidobj))
        if (isempty(avirect))
            F = getframe(gcf);
        else
            F = getframe(gcf, avirect);
        end        
        if (~isempty(aviobj))
            aviobj = addframe(aviobj, F);
        end
        if (~isempty(vidobj))
            writeVideo(vidobj,F);
        end
    end
end
if (nargout > 0)
    varargout{1} = aviobj;
end
end

function datafields(track, ax, inds, TitleOptions, DataPlotOptions, DataAxesOptions)
    fields = {'sbodytheta',track.so.speed_field};
    ttls = {'Body Bend Angle', 'Speed'};
  %  ylb = {'degrees', 'cm/min'};
    mult = [rad2deg(1), 60, 1];
    yl = {[-120 120], [0 mult(2)*max(track.getDerivedQuantity(fields{2}))]};
    for k = 1:2
        hold(ax(k), 'off');
        plot (ax(k), track.dq.eti(inds), mult(k)*track.dq.(fields{k})(inds), DataPlotOptions{:}); hold on
        title(ax(k), ttls{k}, TitleOptions{:});
        set(ax(k), 'YLim', yl{k}, DataAxesOptions{:});   
        %get(ax(k), 'Xlim')
   %     ylabel(ax(k), ylb{k});
    end
  %  set (ax(1), 'YAxisLocation', 'Right');
    if (false)
        ylabel (ax(1), '$\leftarrow$ leftward ; rightward $\rightarrow$', TitleOptions{:}, 'Interpreter', 'Latex');
        ylabel(ax(2), '$\leftarrow$ slower ; faster $\rightarrow$', TitleOptions{:}, 'Interpreter', 'Latex');
    end
end

function handles = updateCenter(handles, ax, track, cind, start, stop, DataAxesOptions)
    if ~isempty(handles)
        for j = 1:length(handles)
            delete(handles(j));
        end
    end
    fields = {'sbodytheta',track.so.speed_field,'vel_dp'};
    mult = [rad2deg(1), 60, 1];
    for k = 1:2
        ih = ishold(ax(k));
        
        hold (ax(k), 'on');
        handles(k) = plot (ax(k), track.dq.eti(cind), mult(k)*track.dq.(fields{k})(cind), 'c.','MarkerSize',30);
        xlim(ax(k), [min(track.dq.eti(start)) max(track.dq.eti(stop))]);
        if (~ih)
            hold (ax(k), 'off');
        end
        if (~isempty(DataAxesOptions))
            set(ax(k), DataAxesOptions{:});
        end
    end
end   

function locPlot (track, inds, ax, varargin)
    field= 'sloc';
    AxesOptions = {};
    imData = [];
    TitleOptions = {};
    varargin = assignApplicable(varargin);
    hold(ax, 'off');
    if (~isempty(imData))
        track.plotPath(field, 'b-', 'inds', inds, 'Axes', ax, varargin{:});
        axis(ax,'equal');
        xl = get(ax, 'XLim');
        yl = get(ax, 'YLim');
       
        u = get(ax, 'Units');
        set(ax,'Units','Pixels');
        pos = get(ax, 'position');
        pos = ceil(pos);
        xinds = find(imData.x <= xl(2) & imData.x >= xl(1));
        yinds = find(imData.y <= yl(2) & imData.y >= yl(1));
        
        if (pos(3) < length(xinds))
            xinds = interp1(xinds, (0:(pos(3)-1)) * length(xinds)/pos(3) + 1, 'nearest');
        end
        if (pos(4) < length(yinds))
            yinds = interp1(yinds, (0:(pos(4)-1)) * length(yinds)/pos(4) + 1, 'nearest');
        end
        set(ax, 'Units', u);
        min(imData.x(xinds))
        max(imData.x(xinds))
        pcolor(ax, imData.x(xinds), imData.y(yinds), imData.im(yinds,xinds)); shading(ax, 'flat'); colormap(ax, 'gray');
        hold(ax, 'on');
    end
    track.plotPath(field, 'b-', 'inds', inds, 'Axes', ax, varargin{:});
    if (exist('xl', 'var'))
        set(ax, 'XLim', xl, 'YLim', yl);
    end
    axis(ax, 'equal');
    title (ax, 'Position', TitleOptions{:});
    hold (ax, 'off');
   
    set(ax,'YDir', 'reverse', AxesOptions{:});
end

function h = locUpdate (track, ind, ax, h, varargin)
    field= 'sloc';
    AxesOptions = {};
    varargin = assignApplicable(varargin);
    existsAndDefault('h', []);
    if (~isempty(h))
        delete(h);
    end
    xl = get(ax, 'XLim');
    yl = get(ax, 'YLim');
    ih = ishold(ax);
    hold(ax, 'on');
   
    loc = track.getDerivedQuantity(field, [], 'inds', ind);
    h = plot(ax, loc(1), loc(2), 'c.', 'MarkerSize', 30);
    if (~ih)
        hold(ax, 'off');
    end
    set(ax, 'XLim', xl, 'YLim', yl, AxesOptions{:});
end



