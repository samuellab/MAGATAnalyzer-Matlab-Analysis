function createGradientsForPhototaxis(cbsize, checkerfn, platefn)
% creates a set of gradients that run along the x-direction of the
% checkerboard
% function createGradientsForPhototaxis(cbsize, checkerfn, platefn)
%
% cbsize is the size of a checker in the output checkerboard image (in cm)
% defaults to 3 cm;
% if checkerfn, platefn, are not passed in, they are selected in a dialog

existsAndDefault('cbsize', 3);
if (~existsAndDefault('checkerfn', []))
    [fn, chbasedir] = uigetfile('*.jpg;*.bmp;*.tiff;', 'Select visual image of checkerboard', fullfile('\\labnas1', 'Share', 'Phototaxis', 'calibrations',''));
    checkerfn = fullfile(chbasedir, fn);
end
if (~existsAndDefault('platefn', []))
    [fn, basedir] = uigetfile('*.jpg;*.bmp;*.tiff;', 'Select visual image of plate', chbasedir);
    platefn = fullfile(basedir, fn);
end

chim = imread(checkerfn);
pim = imread(platefn);

cc = CameraCalibration(chim);

imagesc(pim); colormap (gray(256)); caxis([0 128]);
while(1)
    title ('please click on low edge of plate & hit enter');
    [x,y] = getpts;
    try
        x1 = cc.realPtsFromCamPts([x(1); y(1)]);
    catch
        continue;
    end
    break;
end
lowx = x;
while(1)
    title ('please click on high edge of plate & hit enter');
    [x,y] = getpts;
    try
        x2 = cc.realPtsFromCamPts([x(1); y(1)]);
    catch
        continue;
    end
    break;
end
if (abs(x - lowx) < 0.9 * size(pim,1))
    while(1)
        title ('please click on top and bottom of plate & hit enter');
        [x,y] = getpts;
        if (length(y) < 2)
            continue;
        end
        try
            y = cc.realPtsFromCamPts([x y]');
        catch
            continue;
        end
    end
else
    y = cc.realPtsFromCamPts([mean([x lowx]) mean([x lowx]); 1, size(pim,1)]);
end

title ('thank you - calculating gradients, this may take a moment');
pause(0.1);
newdir = fullfile(chbasedir, 'gradients','');
if (~isdir(newdir))
    mkdir(newdir);
end

style = {'linear', 'exponential', 'power'};
name = {'cam target linear.tiff', 'cam target exponential.tiff', 'cam target square.tiff'};
for j = 1:length(style)
    gi = generateGradientImageForPhototaxis(x1(1), x2(1), cc, 0, style{j}, 2);
    imwrite(gi, fullfile(newdir, name{j}), 'tiff');
end

 gi = generateCheckerImageForPhototaxis(x1(1), x2(1), min(y(2,:)), max(y(2,:)), cc, cbsize, 0);
 imwrite(gi, fullfile(newdir, 'cam target checker.tiff'), 'tiff');

