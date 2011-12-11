function im2 = expandImage (im, scale)
%function im2 = expandImage (im, scale)
%
%takes an image im and stretches it by a factor scale 
%then returns the result in im2

x = (0:(size(im,2) - 1)) * scale;
y = (0:(size(im,1) - 1)) * scale;

x2 = 0:x(end);
y2 = 0:y(end);

[x,y] = meshgrid(x,y);
[x2,y2] = meshgrid(x2,y2);

im2 = interp2(x,y,im,x2,y2,'*linear');


