function shortPath = getshortpath(longPath)
%function shortPath = getshortpath(longPath)
%on a windows PC only, returns the 8 character directory structure leading
%up to the file.
%

%adapted from code posted in MATLAB help
%http://www.mathworks.com/matlabcentral/answers/93932-how-can-i-get-the-short-path-for-a-windows-long-path-using-matlab-7-8-r2009a

if (~ispc)
    shortPath = longPath;
    return;
end

ff = {};
while(~isdir(longPath))
    [longPath, f, ext] = fileparts(longPath);
    ff = [ff [f ext]];
    if (isempty(longPath))
        shortPath = longPath;
        return;
    end
end

fs = actxserver('Scripting.FileSystemObject');

shortPath = fullfile(fs.GetFolder(longPath).ShortPath, ff{end:-1:1});

fs.delete;
