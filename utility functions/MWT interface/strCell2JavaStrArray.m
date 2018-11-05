function sa = strCell2JavaStrArray (sc)
%function sa = strCell2JavaStrArray (sc)
%
%converts a cell of strings into a java array of strings

tf = cellfun(@(s) ischar(s), sc);
sc = sc(tf);
sa = javaArray('java.lang.String', length(sc));
for j = 1:length(sc)
    sa(j) = java.lang.String(sc{j});
end