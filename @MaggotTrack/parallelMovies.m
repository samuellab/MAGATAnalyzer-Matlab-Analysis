function varargout = parallelMovies(track, varargin)
%function parallelMovies(track, varargin)
%
%this function was used to make videos for a presentation, but the videos
%are, in fact, confusing
delayTime = 0.05;
ptbuffer = 3;
for j = 1:length(track)
    track(j).expt.openDataFile;
    fid(j) = track(j).expt.fid;
end
AxesList = [];
npts = [];
TitleOptions = {'Color', 'w', 'FontSize', 14};
AxesOptions = {'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'Box', 'on', 'LineWidth', 2, 'Layer', 'top',...
    'XTick', [], 'YTick', []};
aviobj = [];
avirect = [];
varargin = assignApplicable(varargin);
if (isempty(AxesList))
    order = [1 3 4 2];
    for j = 1:4
        AxesList(j) = subplot(2,2,order(j)); %#ok<AGROW>
    end
end

indset = {[], [], [], []};
tnum = {[], [], [], []};
tind = {[], [], [], []};
offset = 0;
for j = 1:length(track)
    runs = track(j).isrun;
    reos = false(size(runs));
    acchs = false(size(runs));
    rejhs = false(size(runs));

    reos([track(j).reorientation([track(j).reorientation.numHS] > 0).inds]) = true;
    acchs([track(j).headSwing([track(j).headSwing.valid] & [track(j).headSwing.accepted]).inds]) = true;
    rejhs([track(j).headSwing([track(j).headSwing.valid] & ~([track(j).headSwing.accepted])).inds]) = true;

    runs = imerode(runs, ones([1 2*ptbuffer + 1]));
    reos = imdilate(reos, ones([1 2*ptbuffer + 1]));
    nh = [zeros([1 2*ptbuffer]) ones([1 2*ptbuffer+1])];
    acchs = imdilate(acchs, nh);
    rejhs = imdilate(rejhs, ones([1 2*ptbuffer + 1]));

    indmap = track(j).getDerivedQuantity('mapInterpedToPts');

    tind{1} = [tind{1} find(runs)];
    tind{2} = [tind{2} find(acchs)];
    tind{3} = [tind{3} find(rejhs)];
    tind{4} = [tind{4} find(reos)];
    
    
    runinds = offset + (indmap(runs));
    reoinds= offset + (indmap(reos));
    accinds = offset + (indmap(acchs));
    rejinds = offset + (indmap(rejhs));
    indset{1} = [indset{1} runinds];
    indset{2} = [indset{2} accinds];
    indset{3} = [indset{3} rejinds];
    indset{4} = [indset{4} reoinds];
    tnum{1} = [tnum{1} repmat(j, size(runinds))];
    tnum{2} = [tnum{2} repmat(j, size(accinds))];
    tnum{3} = [tnum{3} repmat(j, size(rejinds))];
    tnum{4} = [tnum{4} repmat(j, size(reoinds))];
    offset = offset + length(track(j).pt);
end
if isempty(npts)
    npts = max([length(indset{1}), length(indset{2}), length(indset{3}), length(indset{4})]);
end

ttls = {'Runs', 'Accepted Head Sweeps', 'Rejected Head Sweeps', 'Reorientations'};
pt = [track.pt];
for j = 1:npts
    for k = 1:4
        ind = mod(j-1, length(indset{k})) + 1;
        %ind = indset{k}(ind);
        pt(indset{k}(ind)).drawTrackImage([],'Axes', AxesList(k), 'fid', fid(tnum{k}(ind)), varargin{:});
        shading(AxesList(k), 'interp')
        
        pn = {'XLim', 'YLim'};
        pv = get(AxesList(k), pn);
        pinds = tind{k}(ind) + (-50:0);
        pinds = pinds(pinds > 1);
        
        hold (AxesList(k), 'on');
        track(tnum{k}(ind)).plotPath('sloc', 'b.-','Axes',AxesList(k),'inds', pinds);
        
        title(AxesList(k), ttls{k}, TitleOptions{:});
        set(AxesList(k), pn, pv);
        set(AxesList(k), AxesOptions{:});
        hold(AxesList(k), 'off');
        
    end
    if ~isempty(aviobj)
        if (isempty(avirect))
            F = getframe(gcf);
        else
            F = getframe(gcf, avirect);
        end        
        aviobj = addframe(aviobj, F);
    end
    pause(delayTime);
end

if (nargout > 0)
    varargout{1} = aviobj;
end