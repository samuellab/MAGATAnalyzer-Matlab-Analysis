function addtime(expt, timfname)
%loads timing information to the file and adds it to all tracks
%function addtime(expt, timfname)
%experiment.addtime(timfname)
%
%outputs: none
%inputs:
%EXPT: a member of the Experiment class
%TIMFNAME: the filename of a .tim file
%
%each .tim file has 3 columns; only the first is used here
%column 1:  elapsed time in ms.  If the time wraps to 0, a jump of 65536 ms
%is assumed to have occurred (this is what mightex cameras do)
%column 2: number of triggers sent to the camera by this frame (ignored)
%column 3: time in seconds as recorded by the system at the time the pic
%was taken (ignored)
%
%alternately, a .mdat file with the header Mightex
[ps, fn] = fileparts(expt.fname);
existsAndDefault('timfname', fullfile(ps, [fn '.mdat']));
try 
    [pathstr,name,ext] = fileparts(timfname);
catch me
    disp(me.getReport);
    return;
end

if (strcmp(timfname, 'doNOTloadTIME'))
    return;
end
t = [];
ismdat = false;
if (strcmpi (ext, '.mdat'))
    try
        tempstruct = importdata(timfname);
        col = find(strcmpi(tempstruct.colheaders, 'Mightex_TimeStamp'), 1);
        if (~isempty(col))
            t = tempstruct.data(:,col);
            ismdat = true;
        else
            col = find(strcmpi(tempstruct.colheaders, 'bufnum'), 1);
            col2 = find(strcmpi(tempstruct.colheaders, 'bufnum_time'), 1);
            tpf = diff(tempstruct.data(:,col2))./diff(tempstruct.data(:,col));
            tpm = median(tpf(isfinite(tpf)));
            tpm = mean(tpf(isfinite(tpf) & tpf > 0.8 * tpm & tpf < 1.2 * tpm));
            t = (tempstruct.data(:,col)-tempstruct.data(1,col)) * tpm;
            ismdat = true;
        end
    catch %#ok<CTCH>
        timfname = fullfile(pathstr, [name '.tim']);
        ismdat = false;
    end
    if (isempty(t) || ~exist('col', 'var') || isempty(col)) 
         timfname = fullfile(pathstr, [name '.tim']);
         t = [];
         ismdat = false;
    end
end



if (isempty(t))
    try 
        data = load(timfname);
        if (isempty(data))
            disp (['failed to load timing info from ' timfname]);
            return;
        end
        t = data(:,1);
    catch me
        disp(me.getReport());
        disp ('failed to load timing info -- use add time yourself before running any functions');
        return
    end
end

expt.timfname = timfname;


ind = find(diff(t) < 0, 1, 'first');
while (~isempty(ind))
    t((ind+1):end) = t((ind+1):end) + 65536;
    ind = find(diff(t) < 0, 1, 'first');
end

inds = find(diff(t) < 1);
if ~isempty(inds)
    t(inds+1) = t(inds+1)+1;
end

expt.elapsedTime = (t - t(1)) / 1000;
dt = round(median(diff(t)))/1000; %take interpolation time to nearest ms
expt.dr.interpTime = dt;
[expt.track.dr] = deal(expt.dr);
indx = (1:length(t)) - 1;
for j = 1:length(expt.track)
    expt.track(j).addTime (indx, expt.elapsedTime);
end

if (exist('data', 'var') && size(data,2) > 4)
    et2 = (data(:,5) - data(2,5))/1000;
    [et2,I] = unique(et2);
    et1 = expt.elapsedTime;
    p = polyfit(et2(100:end), et1(I(100:end)), 1);
    et2 = p(1)*et2 + p(2);
    gq = GlobalQuantity();
    gq.xField = 'eti';
    gq.fieldname = 'lightTarget';
    gq.xData = et2';
    gq.yData = data(I,4)';
    expt.addGlobalQuantity(gq);
    dll = deriv(gq.yData, 10) ./ deriv(gq.xData, 10);
    gq.yData = dll;
    gq.fieldname = 'dlightTarget';
    expt.addGlobalQuantity(gq);
    if (size(data,2) > 5)
        gq.fieldname = 'projOutput';
        gq.yData = data(I,6)';
        expt.addGlobalQuantity(gq);
    end
end

if (ismdat)
    addMDatFields(expt, tempstruct, expt.elapsedTime);
end

function addMDatFields(expt, datastruct, timestampbyrow)

header = {'CamTargetLightLevel', 'ProjOutputValue'};
gqname = {'lightTarget', 'projOutput'};
derivPts = [4, 0];
for j = 1:length(header)
    col = find(strcmpi(datastruct.colheaders, header{j}), 1);
    if (~isempty(col))
        ydat = datastruct.data(:,col);
        xdat = timestampbyrow(isfinite(ydat));
        ydat = ydat(isfinite(ydat));
        gq = GlobalQuantity;
        gq.xField = 'eti';
        gq.xData = xdat;
        gq.yData = ydat;
        gq.fieldname = gqname{j};
        gq.derivationMethod = @GlobalQuantity.oneDinterpolation;
        expt.addGlobalQuantity(gq);
        if (derivPts(j) > 0)
            gq.yData = deriv(ydat, derivPts(j))./deriv(xdat, derivPts(j));
            gq.fieldname = ['d' gqname{j}];
            expt.addGlobalQuantity(gq);
        end
    end
end

addMetaDataFields(expt);
    