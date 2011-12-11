function tx = thetaAxis(binsize, type, outputtype)
%function tx = thetaAxis(binsize, type, outputtype)
%
%creates a vector with evenly spaced bin centers from -180 to 180 degrees
%(or -pi to pi radians)
%type is either 'd'(egrees) or 'r'(adians)
%outputtype is either 'd' or 'r'
%both are optional and default to degrees
%examples
%thetaaxis(90, 'degrees', 'degrees') gives [-135 -45 45 135]
%thetaaxis(90, 'degrees', 'radians') gives [-3*pi/4 -pi/4 pi/4 3*pi/4]
%thetaaxis(pi/2, 'radians') gives [-135 -45 45 135]
%
%if the binsize is not a divisor of 180 degrees, the last bin will be
%a different size than all other bins
existsAndDefault('type', 'd');
existsAndDefault('outputtype', 'd');

if (lower(type(1)) == 'd')
    low = -180;
    high = 180;
else
    low = -pi;
    high = pi;
end

tx = (low + binsize/2):binsize:high;

if (lower(type(1)) == 'd' && lower(outputtype(1)) == 'r')
    tx = deg2rad(tx);
end
if (lower(outputtype(1)) == 'd' && lower(type(1)) == 'r')
    tx = rad2deg(tx);
end