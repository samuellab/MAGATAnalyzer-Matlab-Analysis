function markHTInvalid(track, thresh, varargin)
%function markHTInvalid(track, thresh, varargin)
%
%marks as invalid any point where spineDist/spineLength is > thresh
%spineDist = average distance travelled spine segments relative to the the
%midpoints of the spine in frame n & n-1
%spineLength = length of the spine along the contour

debug = false;
varargin = assignApplicable(varargin);

ihtvalid = track.getDerivedQuantity('spineDist')./track.getDerivedQuantity('spineLength') <= thresh;

if (any(~ihtvalid))
    pt = [track.pt];
    et = [pt.et];
    eti = track.getDerivedQuantity('eti');
    htv = interp1(eti, double(ihtvalid), et) > 0.99;
    htv = imerode(htv, ones([1 3])); %knock out point on either side as well
    htv = num2cell(logical(htv & [pt.htValid])); %if a point was previously invalid, it's still invalid
    
    [pt.htValid] = deal(htv{:});
    inds = (find([pt.htValid], 1, 'first')):(find([pt.htValid], 1, 'last'));
    
    track.pt = pt(inds);
    track.npts = length(inds);
    if (debug)
        figure(1); clf();
        for billybob = 1:2
            track.playMovie('axisSize', 1.25 * max(track.getDerivedQuantity('spineLength')), 'delayTime', .25, 'inds', find(~[pt.htValid]));
        end
    end
end