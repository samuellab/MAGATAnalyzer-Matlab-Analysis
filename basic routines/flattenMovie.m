function data = flattenMovie (mov)
%function data = flattenMovie (mov)
%changes the movie mov into a hxwxNframes stack of data scaled from 0 to 1
%first converts each colormap to gray using rgb2gray, then maps the image
%to the colormap
%mov should have the fields cdata [hxw uint8] and colormap [256x3 double]
sz = size(mov(1).cdata);
data = zeros(sz(1),sz(2),length(mov),'uint8');

for j = 1:length(mov)
    %convert the colormap to grayscale
    gm = rgb2gray(mov(j).colormap);
    %convert the colormap to uint8
    gm = uint8 (255*gm(:,1)/max(gm(:,1)));
    %the values in the image range from 0 to 255; we need them to go from 1
    %to 256 so we add 1, but first convert to a 16 bit integer to make room
    data(:,:,j) = gm(uint16(mov(j).cdata)+1);
end
