function addBoundary(expt, boundaryDirImage, boundaryDistImage, lightImage)
%adds boundary information to a checkerboard experiment
%im = boundaryImage is a 2D array with the following properties
%
%*im(y,x) = theta (-pi <= theta <= pi) iff x,y is a point on the boundary of
% two squares, and not on the corner (intersection of two boundaries) 
% in this case, theta is the direction from the light square to the dark
% square
%
%*im(y,x) > 2*pi iff x,y is on the corner of two squares
%
%*im(y,x) < -2*pi iff x,y is in the interior of a square
%
%im = lightingImage is an array with the property
%im(y,x) > 0 if x,y is in the light, and im(y,x) < 0 if x,y is in the dark
%
%boundaryImage and lightingImage must be the size of the original camera
%image
gq = GlobalQuantity();

[x,y] = meshgrid(1:size(boundaryDirImage,2),1:size(boundaryDirImage,1));

%x = reshape(x,1,[]);
%y = reshape(y,1,[]);

expt.calculateDerivedQuantity({'sloc', 'theta'});

gq.xField = 'sloc';
gq.xData(:,:,1) = uint16(x);
gq.xData(:,:,2) = uint16(y);
gq.derivationMethod = @GlobalQuantity.twoDinterpolation;

corners = false(size(boundaryDirImage));
inds = sub2ind (size(corners), round(expt.structvar(:,2)), round(expt.structvar(:,1)));
corners(inds) = true;
corners = imdilate(corners, ones(round(expt.boundaryThickness/2)));

gq.fieldname = 'onboundary';
gq.yData = logical(abs(boundaryDistImage) <= expt.boundaryThickness/2 & ~corners);
expt.addGlobalQuantity(gq);

gq.fieldname = 'ininterior';
gq.yData = logical(abs(boundaryDistImage) > expt.boundaryThickness/2 & ~corners);
expt.addGlobalQuantity(gq);

gq.fieldname = 'boundarytheta';
gq.yData = boundaryDirImage;
expt.addGlobalQuantity(gq);

gq.fieldname = 'inlight';
gq.yData = lightImage;
expt.addGlobalQuantity(gq);

gq.fieldname = 'boundaryDistance';
gq.yData = boundaryDistImage;
%gq.yData(corners) = - boundaryDistImage(corners);
expt.addGlobalQuantity(gq);

gq.xField = 'shead';
gq.fieldname = 'headBoundaryDistance';
expt.addGlobalQuantity(gq);

gq.fieldname = 'headInLight';
gq.yData = lightImage;
expt.addGlobalQuantity(gq);

gq.fieldname = 'headOnBoundary';
gq.yData = logical(abs(boundaryDistImage) <= expt.boundaryThickness/2 & ~corners);
expt.addGlobalQuantity(gq);


%because these images are generally so large (btw 5 & 40 MB depending on
%whether double, int, etc.) we delete all gq fields now, rather than save
%them

expt.globalQuantity = [];

end

