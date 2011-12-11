function [beset, bename, infoset] = createBatchFiles(info, processingParams)
%function beset = createBatchFiles(info, processingParams)
%

[be, pp] = defaultExtractionProcessingParams;
existsAndDefault('processingParams', pp);

[f2p, oe] = createFileToProcess(info,processingParams);
info = info(~oe);
f2p = f2p(~oe);


[gt,~,I] = unique({info.genotype});
beind = 0;
for j = 1:length(gt)
    inf = info(I == j);
    fp = f2p(I == j);
    [et, ~, K] = unique({inf.experimentType});
    for k = 1:length(et)
        beind = beind + 1;
        beset(beind) = be;    %#ok<*AGROW>
        beset(beind).files_to_process = fp(K == k);
        bename{beind} = [gt{j} '_' et{k} '_' datestr(now, 'yyyy-mm-dd_HH-MM') '.bxx'];
        infoset{beind} = inf(K == k);
    end
   
end

