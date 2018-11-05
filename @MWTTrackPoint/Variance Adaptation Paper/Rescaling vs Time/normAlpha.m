function fitstruct = normAlpha (fitstruct, data, timeType, expnum)
% function fitstruct = normAlpha (fitstruct, data, timeType, expnum)
% normalizes alpha(t) to have mean=1

if (length(fitstruct) > 1)
    for j = 1:length(fitstruct)
        fitstruct(j) = normAlpha (fitstruct(j), data, j);
    end
    return;
end
if (existsAndDefault('expnum', []))
    tinds = data.turnExpnum == expnum;
    rinds = data.runExpnum == expnum;
else
    tinds = true(size(data.turnT));
    rinds = true(size(data.runT));
end

if(strcmpi(timeType, 'ton'))
    if(isfield(data, 'runTon'))
        runT = data.runTon;
        turnT = data.turnTon;
    else
        runT = data.runT;
        turnT = data.turnT;
    end
    nr_ton = histcounts(runT(rinds), binEdgesFromCenters(fitstruct.tx_ton));
    nt_ton = histcounts(turnT(tinds), binEdgesFromCenters(fitstruct.tx_ton));
    wton = (nr_ton + nt_ton)/sum(nr_ton + nt_ton);
    wton = wton(:);
    alpha_ton = fitstruct.alpha_ton';
    norm_factor = sum(alpha_ton.*repmat(wton, [1 size(alpha_ton,2)]),1);
    alpha_ton = alpha_ton ./ repmat(norm_factor, [size(alpha_ton,1) 1]);
    fitstruct.alpha_ton = alpha_ton';
else
    nr_t = histcounts(data.runT(rinds), binEdgesFromCenters(fitstruct.tx));
    nt_t = histcounts(data.turnT(tinds), binEdgesFromCenters(fitstruct.tx));
    wt = (nr_t + nt_t)/sum(nr_t + nt_t);
    wt = wt(:);
    alpha = fitstruct.alpha';
    norm_factor = sum(alpha.*repmat(wt, [1 size(alpha,2)]),1);
    alpha = alpha ./ repmat(norm_factor, [size(alpha,1) 1]);
    fitstruct.alpha = alpha';
end

end