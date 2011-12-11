function clist = colorList()
%function clist = colorList()

c = 'bgrcmyk';
s = {'-',':', '-.', '--'};

for j = 1:length(c)*length(s)
    ci = mod(j-1,length(c)) + 1;
    si = ceil(j/length(c));
    clist{j} = [c(ci) s{si}];
end
    
