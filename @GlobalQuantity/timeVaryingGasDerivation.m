function gasval = timeVaryingGasDerivation (xin, xdata, ydata)
%function gasval = timeVaryingGasDerivation (xin, xdata, ydata)
%
% the gas value at time t, location x is given by the following formulae
% dist = (x - xdata.origin) . (xdata.flowdir)
% gasval = gasval0(t - dist/xdata.flowspeed)
% where gasval at 0 is given by ydata vs. xdata.et
try 
    fnames = fieldnames(xin);
    if (isfield(xin, 'eti'))
        eti = xin.eti;
    else
        eti = xin.(fnames{1});
    end
    if (isfield(xin, 'sloc'))
        loc = xin.sloc;
    else
        loc = xin.(fnames{2});
    end

    flowspeed = xdata.flowspeed;
    flowdir = xdata.flowdir;
    origin = xdata.origin;
    et = xdata.et;

    d = (loc(1,:) - origin(1))*flowdir(1) + (loc(2,:) - origin(2))*flowdir(2);
    toff = d/flowspeed;
    gasval = interp1(et, ydata, eti-toff, 'linear', NaN);
    gasval(~isfinite(gasval)) = interp1(et, ydata, (eti(~isfinite(gasval))-toff(~isfinite(gasval))), 'nearest', 'extrap');
catch me
    disp ('error adding time varying gas quantity.  make sure xdata has fields flowspeed, flowdir, origin, and et; make sure xin has a time field followed by a location field');
    disp (me.getReport);
    gasval = [];
end
