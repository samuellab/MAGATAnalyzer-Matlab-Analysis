function mo = makeMovie(track, mo, varargin)
% plays a movie of the track with annotation
% function playMovie(track, varargin)
% outputs, none
% inputs:
%   TRACK < Track;
%   MO < movieOptions;
%   VARARGIN:
%   enter options as pairs, caps matter
%   options, with defaults
%
%   inds = 1:length(track.pt); -- inds to play
%   iinds = []; interped inds;  if passed, we find inds =
%               gdq(mapinterpedtopts,iinds)
%   startLoc = [];
%   stopLoc = []; if startLoc & stopLoc are both not empty, we run the movie
%       between these two points
%   startTime = [];
%   stopTime = []; if startTime & stopTime are both not empty, we run the
%       movie between these times
%


iinds = [];
inds = 1:length(track.pt);
startLoc = [];
stopLoc = [];
startTime = [];
stopTime = [];
fid = [];
varargin = assignApplicable(varargin);
if (isempty(fid))
    track.expt.openDataFile;
    fid = track.expt.fid;
end

%prepare figure window to show movie
if (isempty(mo.figHandle))
   mo.figHandle = gcf;
end
figure(mo.figHandle);
clf(mo.figHandle);
set(mo.figHandle, 'DefaultTextInterpreter', 'Latex', 'DefaultTextFontSize', mo.fontSize);

if ((mo.makeAvi || mo.addToAvi) && mo.makeImStack)
    disp ('You can''t make both an image stack and an avi -- pick one');
    return;
end

if (mo.makeAvi)
    mo.avi = avifile(mo.aviName, mo.avioptions{:});
    mo.addToAvi = true;
end

if (mo.addToAvi)
    p = get(0, 'ScreenSize');
    p(1) = round(p(3)/2 - mo.aviResolution(1)/2);
    p(2) = round(p(4)/2 - mo.aviResolution(2)/2);
    p(3) = mo.aviResolution(1);
    p(4) = mo.aviResolution(2);
    set(mo.figHandle, 'Position', p);
end
if (mo.makeImStack)
    p = get(0, 'ScreenSize');
    p(1) = round(p(3)/2 - mo.jpegResolution(1)/2);
    p(2) = round(p(4)/2 - mo.jpegResolution(2)/2);
    p(3) = mo.jpegResolution(1);
    p(4) = mo.jpegResolution(2);
    set(mo.figHandle, 'Position', p);
end

set(mo.figHandle, 'Color', mo.backgroundColor);
if (isempty(mo.Axes))
    numAxes = (length(mo.datafields) + 1);
    nh = ceil(sqrt(numAxes));
    nv = ceil(numAxes/nh);
    for j = 1:numAxes;
        mo.Axes(j) = subplot(nh,nv,j, 'Parent', mo.figHandle);
        set(mo.Axes(j), 'Color', mo.dataBackgroundColor, 'XColor', mo.dataAxesColor, 'YColor', mo.dataAxesColor,...
            'Alim', [0 1], 'Layer', 'top', 'DrawMode', 'normal','FontSize',mo.fontSize, 'LineWidth', 3, 'Box', 'on');
    end
end
set(mo.Axes(1), 'Color', mo.imBackgroundColor, 'XColor', mo.imAxesColor, 'YColor', mo.imAxesColor);

tlength = length(track.getDerivedQuantity('eti'));
 pt = [track.pt];
if (~isempty(iinds))
    iinds = iinds(iinds >=1 & iinds <= tlength);
    inds = track.getDerivedQuantity('mapinterpedtopts', false,iinds);
else
    if (~isempty(startLoc) && ~isempty(stopLoc))
        [~,s] = track.nearestPoint (startLoc);
        [~,e] = track.nearestPoint (stopLoc);
        if (s > e)
            inds = e:s;
        else
            inds = s:e;
        end
    end

   
    if (~isempty(startTime) && ~isempty(stopTime))
        s = find([pt.et] >= startTime, 1, 'first');
        e = find([pt.et] <= stopTime, 1, 'last');
        inds = s:e;
    end
    inds = inds(inds >= 1 & inds <= length(pt));
    indslim = track.getDerivedQuantity('mapptstointerped', false, inds([1 end]));
    iinds = indslim(1):indslim(2);
end
mo = drawDataFields(track, mo);

handles = [];

if (mo.labelReorientations)
    isreo = false(size(track.isrun));
    isreo([track.reorientation.inds]) = true;
end
loc = track.getDerivedQuantity(mo.locField);    
   
for j = 1:length(iinds)
    ts1 = tic();
    cla(mo.Axes(1));
    pt(inds(j)).drawTrackImage([],'fid', fid, 'Axes', mo.Axes(1), varargin{:}); 
    if (mo.interpImage)
        shading(mo.Axes(1), 'interp');
    end
    set(mo.Axes(1), 'Clim', mo.imCLim);
    hold (mo.Axes(1), 'on');
    sstart = max(iinds(j) - mo.ptbuffer,1);
    send = min(iinds(j) + mo.ptbuffer, tlength);
    plot(mo.Axes(1), loc(1, sstart:send), loc(2, sstart:send), mo.trackStyle, 'LineWidth', mo.trackWidth);
    plot(mo.Axes(1), loc(1, iinds(j)), loc(2, iinds(j)), mo.ptMarker, 'MarkerSize', mo.ptMarkerSize);
    
    %track.plotPath(mo.locField, mo.trackStyle, 'inds', sstart:send, 'Axes', mo.Axes(1), 'LineWidth', mo.trackWidth);
    %track.plotPath(mo.locField, mo.ptMarker, 'inds', iinds(j), 'Axes', mo.Axes(1), 'MarkerSize', mo.ptMarkerSize);
    
    axis (mo.Axes(1), [loc(1,iinds(j)) + [-mo.imAxisSize/2 mo.imAxisSize/2], loc(2,iinds(j)) + [-mo.imAxisSize/2 mo.imAxisSize/2]]);
    axis (mo.Axes(1), 'equal'); 
    axis (mo.Axes(1), [loc(1,iinds(j)) + [-mo.imAxisSize/2 mo.imAxisSize/2], loc(2,iinds(j)) + [-mo.imAxisSize/2 mo.imAxisSize/2]]);
   
    hold (mo.Axes(1), 'off');
    set(mo.Axes(1), 'XTick', [], 'YTick', []);
    mo = track.makeMovieTrackSpecific(mo, pt(inds(j)), iinds(j));
    
    xl = get(mo.Axes(1), 'XLim');   
    yl = get(mo.Axes(1), 'YLim');
    if (mo.labelRuns && ~isempty(track.run) && track.isrun(iinds(j)))
        text(xl(1), yl(1), 'Run', 'Color', mo.runColor, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom',...
            'Parent', mo.Axes(1));      
    end
    if (mo.labelReorientations && ~isempty(track.reorientation) && isreo(iinds(j)))
        text(xl(2), yl(1), 'Reorientation', 'Color', mo.reorientationColor, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom',...
            'Parent', mo.Axes(1));      
    end
    
    handles = updateCenter(track, iinds(j), mo, handles);
    for k = 2:numAxes;
        set(mo.Axes(k), 'Color', mo.dataBackgroundColor, 'XColor', mo.dataAxesColor, 'YColor', mo.dataAxesColor,...
        'Alim', [0 1], 'Layer', 'top', 'DrawMode', 'normal','FontSize',mo.fontSize, 'LineWidth', 1, 'Box', 'on');
    end
    timeleft = mo.delayTime - toc(ts1);
    if (mo.addToAvi)
        mo.avi = addframe(mo.avi, getframe(mo.figHandle));
        pause(0.01);
    else
        if (mo.makeImStack)
            drawnow('expose');
            figure(mo.figHandle);
            saveas(mo.figHandle, [mo.jpegstub '_' num2str(mo.imstackIndex) '.bmp']);
            mo.imstackIndex = mo.imstackIndex + 1;
        else
            if (timeleft > 0)
                pause(timeleft);
            else
                pause(0.01);
            end
        end
    end
   
end
if (mo.makeAvi && mo.closeAviWhenDone)
    mo.avi = close(mo.avi);
    mo.addToAvi = false;
end

function mo = drawDataFields(track, mo)
if (isempty(mo.datafields))
    return;
end
for j = 1:min(length(mo.datafields), (length(mo.Axes)-1))
    Ax = mo.Axes(j+1);
    xdata = track.getDerivedQuantity('eti');
    ydata = mo.dataOps{j}(track.getDerivedQuantity(mo.datafields{j}));
    cla(Ax);
    plot(Ax, xdata, ydata,  mo.dataLineColor, 'LineWidth', mo.dataLineWidth); 
    hold(Ax, 'all');
    if (isempty(mo.dataYLimits{j}))
        mo.dataYLimits{j} = [min(ydata) max(ydata)];
    end
    ypatchvertex = [mo.dataYLimits{j}(1) mo.dataYLimits{j}(1) mo.dataYLimits{j}(2) mo.dataYLimits{j}(2)]';
   
    if (mo.labelRuns)
        xstart = track.run.getDerivedQuantity('eti','position', 'start');
        xend = track.run.getDerivedQuantity('eti','position', 'end');
        patch([xstart;xend;xend;xstart], repmat(ypatchvertex, 1, length(xstart)), mo.runColor, ...
            'faceAlpha', 0.4, 'edgeAlpha', 0, 'Parent', Ax);
    end
   if (mo.labelReorientations)
        xstart = track.reorientation.getDerivedQuantity('eti','position', 'start');
        xend = track.reorientation.getDerivedQuantity('eti','position', 'end');
        patch([xstart;xend;xend;xstart], repmat(ypatchvertex, 1, length(xstart)), mo.reorientationColor, ...
            'faceAlpha', 0.4, 'edgeAlpha', 0, 'Parent', Ax);
    end
    if (~isempty(mo.dataOverlays) && ~isempty(mo.dataOverlays{j}))
        if (~iscell(mo.dataOverlays{j}))
            mo.dataOverlays{j} = mo.dataOverlays(j);
        end
        h = [];
        for k = 1:length(mo.dataOverlays{j})
            h(k) = plot(Ax, xdata, repmat(mo.dataOps{j}(parseDotField(track, mo.dataOverlays{j}{k})), size(xdata)),...
                'LineWidth', mo.dataLineWidth); %#ok<AGROW>
            if (mo.mirrorOverlays(j))
                plot(Ax, xdata, repmat(-mo.dataOps{j}(parseDotField(track, mo.dataOverlays{j}{k})), size(xdata)),...
                 'LineWidth', mo.dataLineWidth, 'Color', get(h(k), 'Color'));
            end
        end
       % legend(h, mo.overlayLabels{j});
    end
             
        
    title(Ax, mo.dataTitles{j}, 'Color', mo.dataAxesColor);
    if (~isempty(mo.dataYLabels{j}))
        ylabel(Ax, mo.dataYLabels{j}, 'Color', mo.dataAxesColor);
    end
    xlabel(Ax, 'time (s)', 'Color', mo.dataAxesColor);
    set(Ax, 'Color', mo.dataBackgroundColor, 'XColor', mo.dataAxesColor, 'YColor', mo.dataAxesColor,...
            'Alim', [0 1], 'Layer', 'top', 'DrawMode', 'fast');
end

function handles = updateCenter(track, cind, mo, handles)
if ~isempty(handles)
    for j = 1:length(handles)
        delete(handles(j));
    end
end
if (isempty(mo.datafields))
    handles = [];
    return;
end
for j = 1:length(mo.datafields)    

    Ax = mo.Axes(j+1);
    xdata = track.getDerivedQuantity('eti',false,cind);
    ydata = mo.dataOps{j}(track.getDerivedQuantity(mo.datafields{j}, false, cind));
    handles(j) = plot (Ax, xdata, ydata,  mo.dataCurrentMarker, 'MarkerSize', mo.dataCurrentMarkerSize);
    xlim(Ax, xdata + [-mo.dataXRange mo.dataXRange]);
    ylim(Ax, mo.dataYLimits{j});
end
    
 


