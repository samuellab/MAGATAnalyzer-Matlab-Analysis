function str = getReport(track, startInd, endInd, varargin)
% returns in str a text report about the segment of track 
% function str = getReport(track, startInd, endInd, varargin)
%
% output: 
%   STR - a string with the report
% input:
%   TRACK < Track
%   startInd: starting index of segment ([] for 1)
%   endInd: ending index of segment ([] for end of track)
% optional args:
%   'tpreport', [true]/false - whether to generate a report for each
%   run & reorientation

existsAndDefault('startInd', 1);
existsAndDefault('endInd', length(track.getDerivedQuantity('eti')));

if (startInd > endInd)
    temp = startInd;
    startInd = endInd;
    endInd = temp;
end

tpreport = true;
precision = '%.1f';
varargin = assignApplicable(varargin);

time = track.getDerivedQuantity('eti', false, [startInd endInd]);
loc = track.getDerivedQuantity('iloc', false, [startInd endInd]);
displacement = diff(track.getDerivedQuantity('pathLength', false, [startInd endInd]));

locStr = ['Segment starts at time ' num2str(time(1), precision) ' x,y = ' num2str(loc(1,1),precision) ',' num2str(loc(2,1),precision)...
            ' and ends at time ' num2str(time(2),precision) ' x,y = ' num2str(loc(1,2), precision) ',' num2str(loc(2,2),precision)];
        
dx =  diff(loc,[],2);

straightDist = sqrt(sum(dx.^2));
straightAngle = rad2deg(atan2(dx(2), dx(1)));
displacementStr = ['Elapsed time: ' num2str(diff(time), precision) 's Straight line distance = ' num2str(straightDist, precision) ...
        ' angle = ' num2str(straightAngle, precision) '(d) Path length = ' num2str(displacement, precision)];
    
    
if (~isempty(track.run))
    runinds = find([track.run.endInd] > startInd & [track.run.startInd] < endInd);
    if ~isempty(runinds)
        runstr = ['Segment intersects ' num2str(length(runinds)) ' runs'];
    else
        runstr = 'Segment intersects no runs';
    end
else
    runstr = 'Track was not segmented, or has no runs';
end

if (~isempty(track.reorientation))
    reoinds = find([track.reorientation.endInd] > startInd & [track.reorientation.startInd] < endInd);
    if ~isempty(reoinds)
        reostr = ['Segment intersects ' num2str(length(reoinds)) ' reorientations'];
    else
        reostr = 'Segment intersects no reorientations';
    end
else
    reostr = [];
end


str = {locStr, displacementStr, runstr, reostr};
if (tpreport)
    if (~isempty(track.run) && ~isempty(track.reorientation))       
        tp = [num2cell(track.run(runinds)), num2cell(track.reorientation(reoinds))];
        si = zeros(size(tp));
        for j = 1:length(tp)
            si(j) = tp{j}.startInd;
        end
        [~,I] = sort(si);
        for j = I
            rp = tp{j}.getReport('precision', precision);
            str = [str, '----', rp{:}];
        end
    else
        if (~isempty(track.run))
            for j = runinds
                rp =  track.run(j).getReport('precision', precision);
                str = [str, '----', rp{:}];
            end
        end
        if (~isempty(track.reorientation))
            for j = reoinds
                rp =  track.reorientation(j).getReport('precision', precision);
                str = [str, '----', rp{:}];
            end
        end
    end
end