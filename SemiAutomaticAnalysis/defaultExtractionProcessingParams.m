function [be,pp] = defaultExtractionProcessingParams ()
%function [be,pp] = defaultExtractionProcessingParams ()
%
% batchExtractor and processingParams are loaded from
% defaultExtractionParams.mat, in this directory
% if defaultExtractionParams.mat does not exist, they are loaded from 
% defaultBatchFile.bxx
%
thisdir = fileparts(mfilename('fullpath'));%(which('defaultExtractionProcessingParams'));
if (exist(fullfile(thisdir, 'defaultExtractionParams.mat'), 'file'))
    load (fullfile(thisdir, 'defaultExtractionParams.mat'), 'be', 'pp');
else
    be = ReadYaml(fullfile(thisdir, 'defaultBatchFile.bxx'));
    pp = be.files_to_process.processing_params;
    save (fullfile(thisdir, 'defaultExtractionParams.mat'), 'be', 'pp');
end
