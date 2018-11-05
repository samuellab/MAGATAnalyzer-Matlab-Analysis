function pts = bezierFromControlPoints(cpts, npts)
%function pts = bezierFromControlPoints(cpts, npts)
%cpts = 2xnum ctrl points
%npts = number of points to interpolate by
%pts = 2xnpts bezier curve

t = linspace(0,1,npts);
M = bezierMatrix(size(cpts,2)-1, t);
pts = cpts*M;

% 
% x = repmat(cpts(1,end:-1:1), [npts 1])';
% y = repmat(cpts(2,end:-1:1), [npts 1])';
% for d = size(x,1):-1:1
%     for j = 1:(d-1)
%         x(j,:) = t.*x(j,:) + (1-t).*x(j+1,:);
%         y(j,:) = t.*y(j,:) + (1-t).*y(j+1,:);
% 
%     end
% end
% pts = squeeze([x(1,:);y(1,:)]);
% 
% 
% end

