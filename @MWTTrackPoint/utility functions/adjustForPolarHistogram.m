function values = adjustForPolarHistogram(values, binCenters)
%function values = adjustForPolarHistogram(values, binCenters)

lowedge = 1.5*binCenters(1) - 0.5*binCenters(2);
values = mod(values - lowedge, 2*pi) + lowedge;
