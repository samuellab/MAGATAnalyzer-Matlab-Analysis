function lns = lineEdgesInImage(im, varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

sigma = 3;
thresh = 0.25;
eraseDist = 15;
debug = false;
varargin = assignApplicable(varargin);
se = strel('disk', eraseDist,4);
[xd,yd] = imgradient(im, sigma);
[xdl, ydl] = imgradient(im, sigma*2); %lowpassed energy image

e = xdl.^2 + ydl.^2;
e = e/max(e(:));
e = e.*(e > thresh);
j = 0;
while (any(e(:) > 0))
    [~,I] = max(e(:));
    
    [y,x] = ind2sub(size(e),I);
    erasemat = ones(size(e));
    erasemat(y,x) = 0;
 
    [x,y] = straightLineEdgeThruPoint([x y], im ,xd, yd);
    if (~isempty(x) && ~isempty(y))
        j = j+1;
        lns(j).x = x;
        lns(j).y = y;
        lns(j).len = sqrt(diff(x).^2 + diff(y).^2);
        xi = round(linspace(x(1),x(2), ceil(lns(j).len)));
        yi = round(linspace(y(1),y(2), ceil(lns(j).len)));
        valid = (xi >= 1 & xi <= size(erasemat,2) & yi >= 1 & yi <= size(erasemat, 1));
        erasemat(yi(valid), xi(valid)) = 0;
        se2 = strel('line', lns(j).len/2,atan2(diff(y), diff(x)));
        erasemat = imerode(erasemat, se2);
    end
    erasemat = imerode(erasemat, se);
    e = e.*(erasemat);    
    if (~exist('lns', 'var'))
        continue;
    end
    if debug
        figure(1);
        imagesc(im); colormap gray; hold on; axis image
        for k = 1:length(lns)
            plot (lns(k).x, lns(k).y, 'r-', 'LineWidth', 2);
        end
        hold off
        figure(2);
        imagesc(e); colormap jet; axis image;
      %  figure(3);
       % imshow(erasemat);
        pause(0.01);
        
    end
end

