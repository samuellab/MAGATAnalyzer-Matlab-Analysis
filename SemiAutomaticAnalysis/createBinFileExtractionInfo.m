function extractionInfo = createBinFileExtractionInfo(basedir, outputdir, extractionSettings)
%function extractionInfo = createBinFileExtractionInfo(basedir, outputdir, extractionSettings)

es.style = 'marc';
es.imdirname = 'Image Data';
es.extdirname = 'Extracted';
es.data_file_extensions = {'.tim', '.tmp', '.dat'};
es.processingParams = [];
if (nargin == 0)
    extractionInfo = es;
    return;
end

if (exist('extractionSettings', 'var') && isstruct(extractionSettings) && ~isempty(extractionSettings))
    fn = fieldnames(extractionSettings);
    for j = 1:length(fn)
        es.(fn{j}) = extractionSettings.(fn{j});
    end
end

fn = recursiveDirectorySearch(basedir, '*.mmf');

for j = 1:length(fn)
    info(j) = getExtraDataFiles (createOutputFileNames(parseDirectoryAndFileName(fn{j}, 'marc'), basedir, outputdir, es.imdirname, es.extdirname), es.data_file_extensions);
%     pause
end

if (isempty(es.processingParams))
    [~,processingParams] = defaultExtractionProcessingParams;
else
    processingParams = es.processingParams;
end
[beset, bename, infoset] = createBatchFiles(info, processingParams);
for j = 1:length(beset)
    extractionInfo(j).batchExtractor = beset(j);
    extractionInfo(j).batchDestination = bename{j};
    extractionInfo(j).fileInfo = infoset{j};
end
    
    
