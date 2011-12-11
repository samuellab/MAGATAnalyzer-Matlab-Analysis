function [f2p, outputexists] = createFileToProcess(info, processing_params)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

originalFieldName_02138.file_stub = 'file stub';
originalFieldName_02138.output_file = 'output file';
originalFieldName_02138.processing_params = 'processing params';

existsAndDefault('processing_params', []);
pp = processing_params;

 

for j = 1:length(info)
    f2p(j).file_stub = info(j).fstub;
    f2p(j).output_file = info(j).outputBinFile;
    if (~isempty(pp))
        pp.fstub = info(j).fstub;
        pp.extension = info(j).extension(info(j).extension ~= '.');
        pp.headerinfoname = info(j).headerFile;
        pp.logName = info(j).logFile;
    end
    f2p(j).processing_params = pp;

    if (nargout >= 2)
        outputexists(j) = (exist(info(j).outputBinFile, 'file') == 2);
    end
    f2p(j).originalFieldName_02138 = originalFieldName_02138;
end

