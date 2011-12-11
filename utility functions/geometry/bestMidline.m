function [midline, c1, c2, inds] = bestMidline(cpts, endpt1, endpt2)
%function [midline, c1, c2, inds] = bestMidline(cpts, enpt1, endpt2)
%
%splits
cpts = double(cpts);
%constraints = [1:length(cpts);[2:length(cpts) 1]];
%dt = DelaunayTri(double(cpts'), constraints');


if prod(size(endpt1)) == 1
   ind1 = endpt1;
else
   [~,ind1] = min(sum((cpts-repmat(endpt1, 1, size(cpts,2))).^2));
 %  ind1 = dt.nearestNeighbor(endpt1(1), endpt1(2));
end
if prod(size(endpt2)) == 1
   ind2 = endpt2;
else
   [~,ind2] = min(sum((cpts-repmat(endpt2, 1, size(cpts,2))).^2));
   %ind2 = dt.nearestNeighbor(endpt2(1), endpt2(2));
end


for step = [10 5 2 1]
    while(1)
        moi = getmoi(ind1 + (-1:1)*step, ind2, cpts);
        [~,j] = min(moi);
        ind1 = ind1 + j - 2;
        moi = getmoi(ind1, ind2 + (-1:1)*step, cpts);
        [~,k] = min(moi);
        ind2 = ind2 + k - 2;
        if (j == 2 && k == 2)
            break;
        end;
        
    end
end

[mid,c1,c2] = splitOutline(cpts,ind1, ind2);
midline = lowpass1d(mid, length(mid)/20);
inds = [ind1 ind2];

function moi = getmoi(endpt1, endpt2,cpts)
endpt1 = mod(endpt1,length(cpts)) + 1;
endpt2 = mod(endpt2,length(cpts)) + 1;
moi = zeros(length(endpt1), length(endpt2));
for j = 1:length(endpt1)
    for k = 1:length(endpt2)
        [mid, c1, c2] = splitOutline(cpts,endpt1(j), endpt2(k));
        mid = lowpass1d(mid, length(mid)/20);
        moi(j,k) = momentOfInteria(mid, c1, c2);
    end
end
