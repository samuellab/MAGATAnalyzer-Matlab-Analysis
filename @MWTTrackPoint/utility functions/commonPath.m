function p = commonPath (p1, p2)
f1 = java.io.File(p1);
f2 = java.io.File(p2);

p1 = char(f1.getCanonicalPath);
p2 = char(f2.getCanonicalPath);

ll = min(length(p1),length(p2));
ind = find(p1(1:ll) ~= p2(1:ll), 1, 'first');
if (isempty(ind))
    ind = ll;
end
p = p1(1:(ind-1));
if (isempty(strfind('\/:', p(end))))
    p = fileparts(p);
else
    p = fileparts(fullfile(p, ''));
end
