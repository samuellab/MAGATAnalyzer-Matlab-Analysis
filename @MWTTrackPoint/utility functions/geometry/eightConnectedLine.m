function [xout,yout] = eightConnectedLine(x,y)

d = ceil(sqrt((diff(x).^2 + diff(y).^2)));
inds = zeros(1,sum(d));
p = 1;

for j = 1:length(d)
    inds(p:(p+d(j)-1)) = linspace(j,j+(d(j)-1)/d(j),d(j));
    p = p+d(j);
end
xout = round(interp1(x,inds));
yout = round(interp1(y,inds));

valid = [true diff(xout)~=0|diff(yout)~=0];
xout = xout(valid);
yout = yout(valid);