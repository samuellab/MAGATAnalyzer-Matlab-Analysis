function [bwim,thresh] = thresholdToTargetArea (im, targetArea, centerLoc)
%function [bwim,thresh] = thresholdToTargetArea (im, targetArea, centerLoc)
%
%finds a threshold thresh, s.t. the area of the region containing centerLoc
%is approx startThresh

low = double(min(im(:)));
high = double(max(im(:)));
startThresh = (double(high) + double(low))/2;

%existsAndDefault('startThresh', 
%gv = sort(im(:),'descend');
%startThresh = interp1(double(gv), targetArea);
%startThresh
existsAndDefault('centerLoc', (size(im'))/2);
if (size(centerLoc,2) == 1)
    centerLoc = centerLoc';
end
if (im(round(centerLoc(2)), round(centerLoc(1))) < startThresh)
 %   warning('TTA:SECONDREG', 'a brighter off-center region appears in this image');
    startThresh = im(round(centerLoc(2)), round(centerLoc(1)))-1;
end

%existsAndDefault('startThresh', (low + high) / 2);


thresh = double(startThresh);
ncalls = 0;
while (low < (high-1))
    props = regionprops(im >= thresh, 'Centroid', 'Area', 'PixelIdxList');
    %area = sum(im1d >= thresh);
    %area = nnz(im1d >= thresh);
    ind = 1;
    %area = props(1).Area;
    for j = 2:length(props)
        %props(j).Centroid
        %centerLoc
        if (sum((props(j).Centroid - centerLoc).^2) < sum((props(ind).Centroid - centerLoc).^2))
            ind = j;
        end
    end
    area = props(ind).Area;
    if (area == targetArea) 
        break
    end
    if (area < targetArea)
        high = thresh;
    else
        low = thresh;
    end
    thresh = round((low + high) / 2);
    ncalls = ncalls + 1;
end
if (ncalls > 8)
    disp(['called ' num2str(ncalls) 'times']);
end
bwim = false(size(im));
bwim(props(ind).PixelIdxList) = true;
%bwim = im >= thresh;
