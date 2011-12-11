function [dx,validinds] = deriv(x,sigma)
%function [dx,validinds] = deriv(x,sigma)

xx = reshape(x,[],1);
dg = reshape(dgausskernel(sigma),[],1);

dx = reshape(conv2(xx,dg,'same'),size(x));

len = ceil(length(dg)/2);
if (2 * len > length(dx))
    validinds = [];
else
    validinds = (len:length(dx)-len);
end