function [ledLookupTable, ledGlobalQuantity] = createLedTableFromBitFile(expt, bitfilename, addToExpt)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
ledLookupTable = repmat(GlobalLookupTable(),0);
existsAndDefault('addToExpt', true);
 
if (~isprop(expt, 'metadata') || isempty(expt.metadata))
    ds = importdata2(expt.timfname);
else
    ds = expt.metadata;
end

% if (~isfield(ds, 'ledFrame'))
%    
%     return;
% end

if (isfield(expt.metadata, 'ledNBytesOut'))
    try
        firstFrame = find(ds.ledFrame == 0, 1, 'first');
        if (isempty(firstFrame))
            x = ds.bufnum;
            y = ds.ledFrame;
            p = polyfit(y(isfinite(y)),x(isfinite(y)),1);
            firstFrame =  find(x == round(p(2)));
        end
        if (isempty(firstFrame))
            firstFrame = find(ds.ledFrame >= 0, 1, 'first');
            firstFrame = firstFrame - ds.ledFrame(firstFrame);
        end
        ind = find(isfinite(ds.ledNBytesOut) & isfinite(ds.bufnum), 1, 'last');
        bytesperframe = round((ds.ledNBytesOut(ind))/(ds.bufnum(ind) - ds.bufnum(firstFrame)));
        bufaxis = (ds.bufnum - ds.bufnum(firstFrame))*bytesperframe;
        valid = isfinite(bufaxis) & isfinite(expt.elapsedTime);
        bufaxis = bufaxis(valid);
        et = expt.elapsedTime(valid);
        
        bitaxis = min(bufaxis):max(bufaxis);

        timaxis = interp1(bufaxis, et, bitaxis);
        fid = fopen(fixFileNameWin(bitfilename), 'rb');
        bits = fread(fid, inf, 'uint8=>double');
        fclose(fid);
        sigma = bytesperframe/3;
    catch me
        disp(me.getReport());
        return;
    end
else
    try
        firstFrame = find(ds.ledFrame == 0, 1, 'first');
        if (isempty(firstFrame))
            x = ds.bufnum;
            y = ds.ledFrame;
            p = polyfit(y(isfinite(y)),x(isfinite(y)),1);
            firstFrame =  find(x == round(p(2)));
        end
        if (isempty(firstFrame))
            firstFrame = find(ds.ledFrame >= 0, 1, 'first');
            firstFrame = firstFrame - ds.ledFrame(firstFrame);
        end
        ind = find(isfinite(ds.ledNBitsOut) & isfinite(ds.bufnum), 1, 'last');
        bitsperframe = round((ds.ledNBitsOut(ind))/(ds.bufnum(ind) - ds.bufnum(firstFrame)));
        bufaxis = (ds.bufnum - ds.bufnum(firstFrame))*bitsperframe;
        valid = isfinite(bufaxis) & isfinite(expt.elapsedTime);
        bufaxis = bufaxis(valid);
        et = expt.elapsedTime(valid);
        bitaxis = min(bufaxis):max(bufaxis);
        timaxis = interp1(bufaxis, et, bitaxis);
        fid = fopen(fixFileNameWin(bitfilename), 'rb');
        bits = fread(fid, inf, 'ubit1=>double');
        fclose(fid);
        sigma = bitsperframe/3;
    catch me
        disp(me.getReport());
        return;
    end
end

allbits = zeros(size(timaxis));
inds = (bitaxis >= 0);
allbits(inds) = bits(1:nnz(inds));

ledLookupTable = GlobalLookupTable();
ledLookupTable.xField = 'eti';
ledLookupTable.fieldname = 'ledFlicker';
ledLookupTable.xData = timaxis;
ledLookupTable.yData = allbits;
ledLookupTable.derivationMethod =  @GlobalLookupTable.averageInPrecedingBin;

ledLookupTable(2) = ledLookupTable(1);
ledLookupTable(2).fieldname = 'ledFlickerDeriv';
%xData is already evenly sampled/interpolated
ledLookupTable(2).yData = deriv(ledLookupTable(1).yData, sigma)/median(diff(ledLookupTable(1).xData)); %picking sigma somewhat arbitrarily so that our derivative spans about 1 frame

ledLookupTable(3) = ledLookupTable(1);
ledLookupTable(3).fieldname = 'ledFlickerDiff';
ledLookupTable(3).yData = [0 diff(ledLookupTable(1).yData)]/median(diff(ledLookupTable(1).xData));

if(addToExpt)
    expt.addGlobalQuantity(ledLookupTable);
end

ledGlobalQuantity = ledLookupTable.toGlobalQuantity(expt, addToExpt);
%{
ledGlobalQuantity = GlobalQuantity();
ledGlobalQuantity.xField = 'eti';
ledGlobalQuantity.fieldname = 'ledFlicker';
ledGlobalQuantity.xData = unique(expt.elapsedTime(:)');
ledGlobalQuantity.xData = ledGlobalQuantity.xData(isfinite(ledGlobalQuantity.xData));
ledGlobalQuantity.yData = GlobalLookupTable.averageInPrecedingBin(ledGlobalQuantity.xData, timaxis, allbits);
%}
end

