function str = getReport(tp, varargin)
%generates a text report about a track part
%function str = getReport(tp, varargin)
%
%output: STR (a string, or cell of strings if multiple trackparts)
%input: TP < TrackPart
%optional arguments:
%   meanFields: default {'speed'}: fields over which to display a mean
%   value

if (length(tp) > 1)
    for j = 1:length(tp)
        str{j} = tp(j).getReport(varargin{:});
    end
    return;
end

meanFields = {'speed'};
precision = '%.1f';
varargin = assignApplicable(varargin);

time = tp.track.getDerivedQuantity('eti', false, [tp.startInd tp.endInd]);
loc = tp.track.getDerivedQuantity('iloc', false, [tp.startInd tp.endInd]);
displacement = diff(tp.track.getDerivedQuantity('pathLength', false, [tp.startInd tp.endInd]));

locStr = {['Start: time = ' num2str(time(1), precision) ' x = ' num2str(loc(1,1), precision) ', y = ' num2str(loc(2,1), precision)],...
            ['End: time = ' num2str(time(2), precision) ' x = ' num2str(loc(1,2), precision) ', y = ' num2str(loc(2,2), precision)]};
        
dx =  diff(loc,[],2);

straightDist = sqrt(sum(dx.^2));
straightAngle = rad2deg(atan2(dx(2), dx(1)));
displacementStr = ['Elapsed time: ' num2str(diff(time), precision) ' s. Straight line distance = ' num2str(straightDist, precision) ...
        ' angle = ' num2str(straightAngle, precision) '(d) Path length = ' num2str(displacement, precision)];
    
    
meanstr = {};

for j = 1:length(meanFields)
    meanstr = [meanstr ['mean ' meanFields{j} ': ' num2str(tp.getDerivedQuantity(meanFields{j}, 'position', 'mean'), precision)]]; %#ok<AGROW>
end
    
str = [locStr, {displacementStr}, meanstr];