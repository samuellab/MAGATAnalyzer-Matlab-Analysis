function str = batchExtractorToString (be)

for j = 1:length(be.files_to_process)
    be.files_to_process(j).processing_params.fstub = be.files_to_process(j).file_stub;
    [d,fn] = fileparts(be.files_to_process(j).output_file);
    be.files_to_process(j).processing_params.headerinfoname = fullfile(d, [fn '_header.txt']);
    be.files_to_process(j).processing_params.logName = fullfile(d, [fn '_log.txt']);
    be.files_to_process(j).processing_params.outputname = be.files_to_process(j).output_file;
end

str = WriteYamlToString(be);
str = strrep(str, 'thresholdScaleImage: null', 'thresholdScaleImage: ""');
str = strrep(str, '[]', '~'); 

str = regexprep(str, '[^-]  file stub:', '\n-\n  file stub:');
str = regexprep(str, '-[ \f\t\v]+file stub:', '-\n  file stub:');
