function mhr = meanThetaRelativeTo(run, dirField)
%function mhr = meanThetaRelativeTo(run, dirField)
%
%gets the average heading over the course of a run, relative to dirField
%
% run < RUN
% dirField - string; run.getDerivedQuantity(dirField) yields 4-quadrant
% angles

vel = run.getDerivedQuantity('vel');
theta = run.getDerivedQuantity(dirField);
rvel = zeros(size(vel));
rvel(1,:) = cos(theta).*vel(1,:) + sin(theta).*vel(2,:);
rvel(2,:) = cos(theta).*vel(2,:) - sin(theta).*vel(1,:);

mhr = atan2(sum(rvel(2,:)), sum(rvel(1,:)));