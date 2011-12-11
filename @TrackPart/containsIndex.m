function tf = containsIndex(tp, ind) 
%returns a vector of logical indices indicating whether ind is in tp(j)
%function tf = containsIndex(tp, ind) 
%
% TP < trackPart
% ind < positive integer 
if (isempty(tp))
    tf = [];
    return;
end
if (length(tp) > 1)
    tf = false(size(tp));
    for j = 1:length(tp)
        tf(j) = tp(j).containsIndex(ind);
    end
    return;
end

tf = ind >= tp.startInd && ind <= tp.endInd;
