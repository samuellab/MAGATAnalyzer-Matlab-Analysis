function [commonpath, p1residual, p2residual] = getCommonPath (p1, p2)
%function commonpath = getCommonPath (p1, p2)
%
%p1 = 'c:\foo\bar\foobar\yourmom\bob.rob'
%p2 = 'c:\foo\bar\foobar\mymom\bob.rob'
%commonpath = 'c:\foo\bar\foobar'

if (~any(p1 == '.'))
    p1 = fullfile(p1, '.');
end
if (~any(p2 == '.'))
    p2 = fullfile(p2, '.');
end


f1 = java.io.File(p1);
f2 = java.io.File(p2);

p1 = char(f1.getCanonicalPath);
p2 = char(f2.getCanonicalPath);

c1 = pathToCellOfStrings(p1);
c2 = pathToCellOfStrings(p2);

minl = min(length(c1), length(c2));
cc1 = c1(1:minl);
cc2 = c2(1:minl);

tf = strcmpi(cc1, cc2);
if (all(tf))
    ind = length(tf);
else
    ind = find(~strcmpi(c1, c2), 1, 'first') - 1;
end

commonpath = fullfile(c1{1:(ind)});
if (ind == length(c1))
    p1residual = [];
else
    p1residual = fullfile(c1{(ind+1):end});
end
if (ind == length(c2))
    p2residual = [];
else
    p2residual = fullfile(c2{(ind+1):end});
end
