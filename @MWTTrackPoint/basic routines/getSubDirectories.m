function dirlist = getSubDirectories(fn)
%function dirlist = getSubDirectories(fn)

d = dir(fn);
d = d(3:end); %get rid of . and ..

d = d([d.isdir]);

if (fn(end) ~= '\')
    fn = [fn '\'];
end
for j = 1:length(d)
    dirlist{j} = [fn d(j).name];
end

