function d = distance (tp1, tp2)
%@TrackPoint
%d = tp1.distance(tp2)
%function d = distance (tp1, tp2)
%distance between the two track points

d = sqrt (sum (tp1.minus(tp2).^2));