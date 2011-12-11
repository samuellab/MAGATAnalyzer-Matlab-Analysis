temperaturefn = '\\labnas1\Share\Marc Data\temperature ramp larvae\20010204\CS1\CS1_.tmp';
if (~exist ('fepprobetemp', 'var') || isempty(fepprobetemp))
    data = load(temperaturefn);
    fepprobetemp = data(21,:);
    figure(1); clf(1)
    smoothtemp = lowpass1D(fepprobetemp, 2000);
    smoothtemp2 = lowpass1d(fepprobetemp, 4000);
    plot (1:length(fepprobetemp), fepprobetemp, 1:length(smoothtemp), smoothtemp, 1:length(smoothtemp2), smoothtemp2);
    dtemp = deriv(smoothtemp, 1);
    figure(2);clf(2)
    plot (dtemp)

end

fn = '\\labnas1\Share\Marc Data\temperature ramp larvae\20010204\CS1\CS1_tracks_2.bin';
timfn = '\\labnas1\Share\Marc Data\temperature ramp larvae\20010204\CS1\CS1_.tim';

if (~exist ('expt', 'var') || isempty ('expt'))
    expt = Experiment.fromFile(fn, timfn, true, [], 200);
    expt.executeTrackFunction('fixHTOrientation');
end

if (isempty([expt.track.run]) || exist ('startover', 'var') && startover)
    so = MaggotSegmentOptions();
    expt.so = so;
    if (exist ('startover', 'var') && startover)
        expt.executeTrackFunction('recalculateDerivedQuantities');
    end
    expt.executeTrackFunction('segmentTrack', so);
    startover = false;
end

hstime = [];
for j = 1:length(expt.track)
    hstime = [hstime expt.track(j).dq.eti([expt.track(j).headSwing.startInd])];
end
timex = 0:60:max(expt.elapsedTime);
plot (timex, hist(hstime,timex)./hist(expt.gatherField('eti'), timex));
