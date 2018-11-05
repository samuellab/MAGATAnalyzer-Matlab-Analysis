function lpdata = lowpass1D (data, sigma, varargin)
%function lpdata = lowpass1D (data, sigma, varargin)
%lowpasses data with a gaussian filter of width sigma
%
%optional parameter/value pairs
%'padType', pt
%pt is one of 'zeros', 'ends', 'linear', 'mirror'
%
%if padType is zeros, pad ends with zeros
%if padType is 'ends', pad with the values at the end points
%if padType is 'linear', pad(1:kw) = data(1:kw) -
%data(kw) + data(1)
%if padType is 'mirror', pad(1:kw) = 2*data(1) - data(kw:1)
%where kw is kernel width
%if padType is 'circular', pad(1:kw) = data((end-kw+1):end)
%
%if data is a matrix, lowpass is done along the first dimension
%unless the first dimensions is smaller than the second AND also smaller
%than the maximum of (3,kernel size)
%(implying you probably meant to lowpass along the second dimension)
%
% this is a change from previous implementation, and could cause problems
% for some previous code, especially with large first dimensions that are
% meant to be transposed. EG size of the data is (15,2000) and sigma = 1;
% previously would have been transposed, now is not
padType = 'ends';
varargin = assignApplicable(varargin);
g = gaussKernel(sigma);
kw = floor(length(g)/2);
%old code -- this is a bug!
% if (size(data,1) < size(data,2))
%     data = data';
%     transpose = 1;
% else
%     transpose = 0;
% end
%only transpose if first dimension is smaller than second AND
%first dimension is smaller than convolution kernel size
if (size(data,1) < size(data,2) && size(data,1) < max(3,length(g))) 
    data = data';
    transpose = 1;
else
    transpose = 0;
end
lpdata = zeros(size(data));
for k = 1:size(data,2)
    paddeddata = zeros([size(data,1)+2*kw,1]);
    paddeddata((kw+1):end-kw) = data(:,k);
    switch (padType)
        case ('ends')
            paddeddata(1:kw)=data(1, k);
            paddeddata((end-kw+1):end) = data(end,k);
        case ('linear')
            paddeddata(1:kw)=data(1:kw, k) + data(1,k) - data(kw,k);
            paddeddata((end - kw + 1):end) = data(end - kw + 1:end,k) + data(end,k) - data(end - kw + 1, k);
        case ('circular')
            paddeddata(1:kw)=data((end-kw+1):end, k);
            paddeddata((end - kw + 1):end) = data(1:kw,k);
        case ('mirror')
            paddeddata(1:kw)=2*data(1,k) - data(kw:-1:1, k);
            paddeddata((end - kw + 1):end) = 2*data(end,k) - data(end:-1:(end - kw + 1),k);
    end
    lpdata(:,k) = conv2(paddeddata,g','valid');
end
if (transpose)
    lpdata = lpdata';
end