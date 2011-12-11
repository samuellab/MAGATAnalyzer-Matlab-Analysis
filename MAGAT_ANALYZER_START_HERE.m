% this script sets up directories then runs the MAGAT_ANALYZER_DEMO
% for more details see MAGAT_ANALYZER_DEMO
% this script must be run from the Matlab-Track-Analysis directory

setupDirectories
addpath (genpath(fullfile(pwd, 'MAGATAnalyzer Example Scripts')));
help MAGAT_ANALYZER_DEMO
MAGAT_ANALYZER_DEMO