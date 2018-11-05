function diagIm = diagnosticImage (expt, foregroundIm, varargin)

cc = expt.camcalinfo;
if ~existsAndDefault('foregroundIm', []);
    [dr, n] = fileparts(expt.fname);
    d = dir(fullfile(dr, [n(1:end-1) '*.bmp']));
    
    if (isempty(d))
        warning( ['Could not find foreground image corresponding to: ' expt.fname]);
        diagIm = [];
        return;
    end
    im = double(imread(fullfile(dr, d(1).name)));
else
    if (ischar(foregroundIm))
        im = double(imread(foregroundIm));
    else
        im = foregroundIm;
    end
end
pts = expt.gatherField('iloc');
ipts = round(cc.camPtsFromRealPts(pts));
ipts = ipts(:, ipts(1,:) >= 1 & ipts(1,:) <= size(im,2) & ipts(2,:) >= 1 & ipts(2,:) <= size(im,1));
inds = sub2ind(size(im), ipts(2,:), ipts(1,:));
ispath = zeros(size(im));
ispath(inds) = 1;
ispath = imdilate(ispath, ones(ceil(sqrt(mean(expt.gatherField('iarea'))))));
diagIm = repmat(im, [1 1 3]);
diagIm(:,:,2) = im.*double(ispath);
diagIm(:,:,3) = im.*double(ispath);
        