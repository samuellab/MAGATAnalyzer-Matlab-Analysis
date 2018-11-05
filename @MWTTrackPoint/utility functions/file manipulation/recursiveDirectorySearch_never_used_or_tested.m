function fileList = recursiveDirectorySearch(dirName, pattern, maxDepth)
%function fileList = recursiveDirectorySearch(dirName, pattern, maxDepth)
%
%modified from code posted by gnovice and Peter D on stackoverflow
%http://stackoverflow.com/questions/2652630/how-to-get-all-files-under-a-specific-directory-in-matlab

  if (nargin < 3)
      maxDepth = 5; %by default, do not search farther than 5 directories deep
  end
  if (nargin < 2)
      pattern = '';
  end
  if (maxDepth < 0)
      fileList = {};
      return;
  end
  
  dirData = dir(dirName);      %# Get the data for the current directory
  dirIndex = [dirData.isdir];  %# Find the index for directories
  fileList = {dirData(~dirIndex).name}';  %'# Get a list of the files
  
  if (~isempty(pattern) && ~isempty(fileList)) 
      matchstart = regexp(fileList, pattern); 
      fileList = fileList(~cellfun(@isempty, matchstart));
  end
      
  
  if ~isempty(fileList)
    fileList = cellfun(@(x) fullfile(dirName,x),...  %# Prepend path to files
                       fileList,'UniformOutput',false);
  end
  subDirs = {dirData(dirIndex).name};  %# Get a list of the subdirectories
  validIndex = ~ismember(subDirs,{'.','..'});  %# Find index of subdirectories
                                               %#   that are not '.' or '..'
  for iDir = find(validIndex)                  %# Loop over valid subdirectories
    nextDir = fullfile(dirName,subDirs{iDir});    %# Get the subdirectory path
    fileList = [fileList; recursiveDirectorySearch(nextDir, pattern, maxDepth - 1)];  %# Recursively call getAllFiles
  end

end