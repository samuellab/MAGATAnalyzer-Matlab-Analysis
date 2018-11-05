function c = nthcolor (n)
%function c = nthcolor (n)

clist = colorList;
c = clist{mod(n-1,length(clist))+1};

