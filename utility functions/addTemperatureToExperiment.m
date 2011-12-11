function addTemperatureToExperiment(expt, varargin)
%function addTemperatureToExperiment(expt, varargin)
%
% expt < Experiment
% optional arguments:
%       'tmpfname' - the name of the .tmp file; default is to use the
%         expt.fname with the extension changed to .tmp
%       'tmpchannel' - the channel number of the temperature sensor (0 -
%        15) - default is to show a graph and allow the user to enter

if length(expt) > 1
    for j = 1:length(expt)
        addTemperatureToExperiment(expt(j), varargin{:});
    end
    return;
end

[pathstr,name] = fileparts(expt.fname);
tmpfname = fullfile(pathstr, [name '.tmp']);

tmpchannel = [];

varargin = assignApplicable(varargin);

try 
    data = load(tmpfname);
    if (isempty(data))
        disp (['failed to load temperature data from ' tmpfname]);
        return;
    end
catch me
    disp(me.getReport);
    disp (['failed to load temperature data from ' tmpfname]);
    return;
end

chnum = data(1:3:end, 1);
temp = data(3:3:end, :);
tim = data(2:3:end, :);
tim = tim - min(min(tim));
chname = cellfun(@num2str, num2cell(chnum), 'UniformOutput', false);
if ~isempty(tmpchannel)
    ch = find(chnum == tmpchannel);
else
    ch = [];
end

while isempty(ch)
    clf()
    plot (tim', temp');
    legend(chname);
    tmpchannel = input ('input the channel number of the surface temperature probe');
    ch = find(chnum == tmpchannel);
end

[tim,I] = unique(tim(ch,:));
tim = (tim - tim(1))/1000; %convert to seconds
temp = temp(ch,I);

%rescale timing to match time in experiment (.tim file) -- this is because
%the clock on the temperature controller may not run at exactly the same
%speed as the clock on the camera
% report if the difference is greater than 1%

totaltimetc = tim(end);
totaltimecam = expt.elapsedTime(end);

if (abs(totaltimetc - totaltimecam)/totaltimecam > 0.01)
    disp (['warning: camera time differs from temperature controller time by ' num2str(100*abs(totaltimetc - totaltimecam)/totaltimecam) '%']);
    disp(['camera time = ' num2str(totlatimecam) ' and tc time = ' num2str(totaltimetc)]);
    disp (['.tim file = ' expt.timfname]);
    disp (['.tmp file = ' tmpfname]);
end

t = tim * totaltimecam/totaltimetc;
dt = min(min(diff(tim)), expt.dr.interpTime);
ti = 0:dt:max(t);
tempi = interp1(tim,temp,ti);
temps = lowpass1D(tempi, max(30/dt,3));
dtemp = deriv(temps, max(5/dt,3));

gq = GlobalQuantity();
gq.xData = ti;
gq.xField = 'eti';
gq.derivationMethod = @GlobalQuantity.oneDinterpolation;

gq.fieldname = 'temp';
gq.yData = temps;

expt.addGlobalQuantity(gq);

gq.fieldname = 'dtemp';
gq.yData = dtemp;

expt.addGlobalQuantity(gq);

clf();
subplot(2,1,1); plot (ti, temps); title ('temperature vs. time');
subplot(2,1,2); plot (ti, dtemp); title ('temperature derivative vs. time');

