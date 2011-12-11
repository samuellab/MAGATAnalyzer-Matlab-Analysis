function [dx,validinds] = deriv(x,sigma,varargin)
%function [dx,validinds] = deriv(x,sigma,varargin)
%
%optional args:
%padtype can be 'linear' (default) or 'circular' (connects ends)

padtype = 'linear';
padType = [];
varargin = assignApplicable(varargin);
if (~isempty(padType))
    padtype = padType;
end
if (size(x,1) > size(x,2))
    xx = x';
    t = 1;
else
    xx = x;
    t = 0;
end
dg = reshape(dgausskernel(sigma),1,[]);

padfront = ceil ((length(dg)-1) / 2);
padback = length(dg) - padfront - 1;
if (padfront >= size(xx,2)  || padback >= size(xx,2))
    warning ('deriv:sigtoobig', 'deriv called on data set too short to support sigma');
    ks = floor((size(xx(:,2)) - 1)/2);
    inds = ceil(length(dg)/2 + -ks:ks);
    dg = dg(:,inds);
    padfront = ceil ((length(dg)-1) / 2);
    padback = length(dg) - padfront - 1;
end

    
if (strcmpi(padtype, 'circular'))
    frontpad = xx(:,(end-padfront+1):end);
    backpad = xx(:,1:padback);
else
    frontpad = 2*repmat(xx(:,1), 1, padfront) - xx(:,(padfront+1):-1:2);
    backpad = 2*repmat(xx(:,end), 1, padback) - xx(:,(end-1):-1:(end-padback));
end
%xx = [repmat(xx(:,1), 1, padfront) xx repmat(xx(:,end), 1, padback)];
xx = [frontpad xx backpad];
dx = conv2(1,dg,xx,'valid');

if (t)
    dx = dx';
end
len = ceil(length(dg)/2);
if (2 * len > length(dx))
    validinds = [];
else
    validinds = (len:length(dx)-len);
end