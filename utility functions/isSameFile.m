function tf = isSameFile (path1, path2)
%function tf = isSameFile (path1, path2)

if (~ischar(path1) || ~ischar(path2))
    warning ('can only compare string paths');
    tf = false;
    return;
end

if (~any(path1 == '.'))
    path1 = fullfile(path1, '.');
end
if (~any(path2 == '.'))
    path2 = fullfile(path2, '.');
end

f1 = java.io.File(path1);
f2 = java.io.File(path2);

p1 = f1.getCanonicalFile;
p2 = f2.getCanonicalFile;

tf = ~logical(p1.compareTo(p2));
