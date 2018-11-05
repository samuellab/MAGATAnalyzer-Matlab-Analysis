function cpts = removeSmallLoops(cpts)
%function cpts = removeSmallLoops(cpts)
%
%removes any small loops (where same point appears twice in a contour) from
%a contour; it's ok if first and last point are the same (contour is
%closed)
%
%cpts is 2xNpts
npts = size(cpts,2);
cpmat1 = repmat(cpts(1,:), npts, 1);
cpmat2 = repmat(cpts(2,:), npts, 1);
[j,k] = find(cpmat1 == cpmat1' & cpmat2 == cpmat2' & triu(true(size(cpmat1)),+1));

valid = (j~=1 | k~=npts);
j = j(valid); k = k(valid);

valid = true(1,npts);
for m = 1:length(j)
    if (k(m) - j(m)) < npts/2
        valid((j(m)+1):k(m)) = false;
    else
        valid([(k(m)+1):npts 1:j(m)]) = false;
    end
end

cpts = cpts(:,valid);
