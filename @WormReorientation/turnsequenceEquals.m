function tf = turnsequenceEquals(reo, turnsequence)
% true iff reo.trunsequence == turnsequence; works on vectors of reo
% element by element
% function tf = turnsequenceEquals(reo, turnsequence)
if (isempty(reo))
    tf = [];
    return;
end

if (length(reo) > 1)
    tf = false(size(reo));
    for j = 1:length(reo)
        tf(j) = reo(j).turnsequenceEquals(turnsequence);
    end
    return;
end

tf = (length(reo.turnsequence) == length(turnsequence) && all(reo.turnsequence == turnsequence));
