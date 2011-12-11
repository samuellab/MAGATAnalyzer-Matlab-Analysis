function im = makeWholeFrameImage (expt, ind, varargin)
%creates a mosaic image of all the individual point images for a frame
%function im = makeWholeFrameImage (expt, ind, varargin)
%
%outputs:
%IM: a grayscale image
%inputs:
%EXPT: a member of the experiment class
%IND: the frame to make an image of
%optional parameter/value pairs:
%'imsize', [height width], the size of the output image; default [1944 2592] 
%experiment frame ind

imsize = [1944 2592];

assignApplicable(varargin);

im = zeros(imsize);

for j = 1:length(expt.track)
    ptind = find([expt.track(j).pt.ind] == ind);
    if (~isempty(ptind))
        pt = expt.track(j).pt(ptind);
        if (isempty(pt.imData))
            pt = expt.reloadPoint(pt);
        end
        im(double(pt.imOffset(2)) + (1:size(pt.imData,1)), double(pt.imOffset(1)) + (1:size(pt.imData,2))) = double(pt.imData);
    end
end