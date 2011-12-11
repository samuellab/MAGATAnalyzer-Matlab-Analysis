function illustrateHeadSwing (hs, varargin)
%function illustrateHeadSwing (hs, varargin)

position = [];

varargin = assignApplicable(varargin);


theta = hs.track.getDerivedQuantity('sbodytheta', false, 1:(hs.endInd-1));
startInd = find(abs(theta) < hs.track.so.headswing_stop | sign(theta) ~= hs.sign, 1, 'last');
if isempty(startInd)
    startInd = hs.startInd;
end
[blah,maxInd] = max(abs(theta(startInd:(hs.endInd-1))));
if(isempty(maxInd))
    maxInd = hs.endInd - startInd + 1;
end
sind = hs.track.getDerivedQuantity('mapInterpedToPts', false, startInd - 1 + [1 maxInd]);

if (~isempty(position))
    offset = position - mean([hs.track.pt(sind).loc],2);
else
    offset = [0;0];
end

ih = ishold;
trackinds = (hs.startInd - 15):(hs.endInd + 15);
trackinds = trackinds(trackinds > 0 & trackinds < hs.track.npts);
loc = hs.track.getDerivedQuantity('sloc');

loc = loc + repmat(offset, [1 length(loc)]);
plot (loc(1,trackinds), loc(2,trackinds), 'k.','MarkerSize',3);hold on

drawContourAndHead(hs.track.pt(sind(1)), 'offset', offset, 'contourColor', 'k-', varargin{:}); 
if (hs.accepted)
    color = 'g-';
else
    color = 'r-';
end
plot (loc(1,trackinds(trackinds>hs.endInd)), loc(2,trackinds(trackinds>hs.endInd)), [color(1) '.'],'MarkerSize',4);hold on
drawContourAndHead(hs.track.pt(sind(2)), 'offset', offset, 'contourColor', color, varargin{:});

if(~ih)
    hold off
end