function [pt0,x,y,theta] = lightDirectionFromLines(lines, predominantTheta, minLen)
%function [x,y,theta] = lightDirectionFromLines(lines, predominantTheta, minLen)
%takes the line structure from lineEdgesInImage and produces a set of x,y
%positions and the angle of light at the position
%
%predominant theta is the direction the light is coming from in general,
%and is used to resolve ambiguity without having to locate the pins
%
%pt0 is the point from which the light originates

existsAndDefault('predominantTheta', 0);
existsAndDefault('minLen', 100);

lines = lines([lines.len] > minLen);
x = zeros(size(lines));
y = x;
theta = x;

for j = 1:length(lines)
    th = atan2(diff(lines(j).y), diff(lines(j).x));
    dth = diff(unwrap([predominantTheta; th]));
    if (abs(dth) > pi/2)
        x(j) = lines(j).x(2);
        y(j) = lines(j).y(2);
        theta(j) = th + pi;
    else
        x(j) = lines(j).x(1);
        y(j) = lines(j).y(1);
        theta(j) = th;
    end
end

theta = mod(theta + pi - predominantTheta, 2*pi) - pi + predominantTheta; %recenter around predominantTheta

distsqtolines = @(pt) sum((-sin(theta).*(x-pt(1)) + cos(theta).*(y-pt(2))).^2);

pt0 = [x(1) y(1)];
pt0 = fminsearch(distsqtolines, pt0);
