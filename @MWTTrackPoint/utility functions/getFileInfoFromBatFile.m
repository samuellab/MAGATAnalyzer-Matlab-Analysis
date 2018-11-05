function [directories, fstubs, outnames] = getFileInfoFromBatFile(str)
%function [directories, fstubs, extensions, outnames] = getFileInfoFromBatFile(str)
%
%str contains the entire text of the bat file

fsinds = strfind(str, 'file stub:');
ofinds = strfind(str, 'output file:');

for j = 1:length(fsinds)
    stub = strtrim(sscanf(str(fsinds(j):end), 'file stub:%[^\n]s'));
    i = strfind(stub, '\');
    i = i(end);
    directories{j} = stub(1:i);
    fstubs{j} = stub((i+1):end);
    
    outnames{j} = strtrim(sscanf(str(ofinds(j):end), 'output file:%[^\n]s'));
end