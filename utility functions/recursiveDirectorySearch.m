function fileNames = recursiveDirectorySearch(startingDirectory, fileSpecifier)
%function fileNames = recursiveDirectorySearch(startingDirectory, fileSpecifier)
%
%returns complete paths to all files in startingDirectory or subdirectories
%that satisfy fileSpecifier
%
%example
%startingDirectory = '\\labnas1\share\David\Extracted\Spatial'
%batfiles = recursiveDirectorySearch(startingDirectory, '*.bat')
if (ispc)
    cmd = ['dir /s /b "' fullfile(startingDirectory, '*.mmf') '"'];
    
else
    if (fileSpecifier(1) == '*')
        fileSpecifier = ['\' fileSpecifier];
    end
    cmd = ['find "' startingDirectory '" -name ' fileSpecifier];
%    [s,w] = system(); % this has not been tested yet MHG 7/6/2011
end
[s,w] = system(cmd);
if (s ~= 0)
    warning (['system command ' cmd ' returned an error code']);
end
fileNames = regexp(w, '[\r\n]+', 'split');

fileNames = fileNames(cellfun(@(s) ~isempty(s),fileNames));
%{
p = genpath(startingDirectory);
inds = [0 strfind(p,';')];
for j = 1:(length(inds) - 1)
    directories{j} = p((inds(j)+1):(inds(j+1) - 1));
end

fnnum = 0;
for j = 1:length(directories)
    d = dir([directories{j} '\' fileSpecifier]);
    for k = 1:length(d)
        fnnum = fnnum+1;
        fileNames{fnnum} = [directories{j} '\' d(k).name];
    end
end
%}
% 
% function shortPath = getshortpath(longPath)
% fs = actxserver('Scripting.FileSystemObject');
% shortPath = fs.GetFolder(longPath).ShortPath;
% fs.delete;
% 
% function longPath = getlongpath(shortPath)
% fs = actxserver('Scripting.FileSystemObject');
% if (iscell(shortPath))
%     for j = 1:length(shortPath)
%         longPath{j} = fs.GetFolder(shortPath{j}).Path;
%     end
% else
%     longPath = fs.GetFolder(shortPath).Path;
% end