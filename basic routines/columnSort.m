function [x,y] = columnSort (x,y, ncols, rorder)
%function [x,y] = columnSort (x,y, ncols, rorder)
%
%sorts a set of x,y points into row, column order
%rorder is either 'ascend' in which case the first row is low y values
%or 'descend' in which case the first row has high y values
%(default 'ascend')

existsAndDefault('rorder', 'ascend');
[y,I] = sort(y, rorder);
x = x(I);

nrows = length(x)/ncols;

 for j = 1:nrows
     inds = (j-1)*ncols + (1:ncols);
     [xt,I] = sort(x(inds));
     row(j).x = xt;
     row(j).y = y(inds(I));
 end
 x = [row.x];
 y = [row.y];