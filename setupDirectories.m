function setupDirectories (username)
%function setupDirectories (username)
%adds paths to the matlab path so everything will run correctly
%(@Class directories do not need to be added to the path)
%if a username is given, adds user specific\username and subdirectories
%to the path
addpath (pwd);
addpath (fullfile(pwd, 'basic routines'));
addpath (fullfile(pwd,'example scripts'));
addpath (fullfile(pwd,'useful extra classes'));
if (exist (fullfile (pwd, 'MartAnalysis'),'dir'))
    addpath (fullfile(pwd,'MartaAnalysis'));
end
if (exist (fullfile(pwd, 'AndyAnalysis'), 'dir'))
    addpath (fullfile(pwd, 'AndyAnalysis'));
    addpath (fullfile(pwd, 'AndyAnalysis', 'YMLUtils'));
end
addpath (genpath(fullfile(pwd,'utility functions')));
if (exist(fullfile(pwd, 'guis'), 'dir'))
    addpath (genpath(fullfile(pwd,'guis')));
end
if (exist (fullfile(pwd, 'Simulation Classes'), 'dir'))
    addpath (fullfile(pwd,'Simulation Classes'));
end
addpath([docroot '/techdoc/creating_plots/examples']);
if (exist ('username', 'var') && ~isempty(username))
    addpath (genpath(fullfile(pwd, 'user specific', username, '')));
end
addpath(genpath(fullfile(pwd, 'yamlMatlab')));
addpath(fullfile(pwd, 'SemiAutomaticAnalysis'));
