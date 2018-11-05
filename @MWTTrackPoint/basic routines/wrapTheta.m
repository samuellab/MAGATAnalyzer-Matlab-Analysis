function theta2 = wrapTheta (theta)
%function theta2 = wrapTheta (theta)
%
%adds factors of 2*pi to eliminate large jumps in theta

theta2 = theta;
dt = diff([theta(1) theta]);
%size(dt)
indsup = find(dt < -3/2*pi);
indsdown = find(dt > 3/2*pi);

%plot(1:length(theta),theta,'b-',indsup,theta(indsup),'r.',indsdown,theta(indsdown),'g.');

for j = indsup
    theta2(j:end) = theta2(j:end) + 2*pi;
end

for j = indsdown
    theta2(j:end) = theta2(j:end) - 2*pi;
end

%plot(1:length(theta),theta2,'b-',indsup,theta2(indsup),'r.',indsdown,theta2(indsdown),'g.');


