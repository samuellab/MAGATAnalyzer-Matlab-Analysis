function calculateDerivedQuantity(track, quantityNames, recalculate)
%calculates derived quantity(ies) and stores in track.dq
%function calculateDerivedQuantity(track, quantityNames, recalculate)
%
%inputs: 
%TRACK: a member of the Track class
%QUANTITYNAMES: a list of quantities to calculate; for a list of valid
%   quantities, see Track/validDQName
%RECALCULATE: if true, calculates the quantity again even if it already
%   exists as a field of track.dq (default false)
    if ~iscell(quantityNames)
        quantityNames = {quantityNames};
    end
    if (~exist('recalculate', 'var'))
        recalculate = false;
    end
    for j = 1:length(quantityNames)
        if (isfield(track.dq, quantityNames{j}) && ~recalculate)
            continue;
        end
        switch(quantityNames{j})
            case 'eti'
                calculateInterpedTime(track);
            case 'etiFromTrackStart'
                track.calculateDerivedQuantity('eti');
                track.dq.etiFromTrackStart = track.dq.eti - track.dq.eti(1);
            case 'iloc'
                track.calculateDerivedQuantity('eti', false);
                calculateInterpedLocation(track);
            case 'sloc'
                track.calculateDerivedQuantity('iloc', false);
                calculateSmoothedLocation(track);
            case {'xloc', 'yloc'}
                track.calculateDerivedQuantity('sloc', recalculate);
                track.dq.xloc = track.dq.sloc(1,:);
                track.dq.yloc = track.dq.sloc(2,:);
            case {'vel', 'speed', 'vnorm', 'theta', 'nspeed', 'nvel'}
                track.calculateDerivedQuantity({'iloc','sloc'}, false);
                caclulateVelocity(track);
            case {'adjspeed'}
                track.calculateDerivedQuantity('speed', false);
                calculateAdjSpeed(track);
            case {'speed_diff_local'}
                track.calculateDerivedQuantity('speed');
                calculateSpeedDiff(track);
            case {'deltatheta', 'ddtheta', 'acc', 'curv'}
                track.calculateDerivedQuantity({'vel', 'theta'}, false);
                calculateAcceleration(track);
            case {'lrdtheta'}
                calculatelrdtheta(track);                
            case {'pathLength'}
                track.calculateDerivedQuantity({'iloc','sloc'}, false);
                calculatePathLength(track);
            case {'displacement'}
                track.calculateDerivedQuantity({'iloc','sloc'}, false);
                calculatePathLength(track);
            case {'icov', 'covRatio', 'covTheta', 'covMinor', 'covMajor'}
                track.calculateDerivedQuantity('eti', false);
                calculateCovariance(track);
            case {'scov', 'scovRatio', 'scovTheta', 'scovMinor', 'scovMajor'}
                track.calculateDerivedQuantity('icov', false);
                calculateSmoothedCovariance(track);
            case {'dcovRatio'}
                calculateDerivativeOfCovarianceRatio(track);
            case {'iarea', 'sarea'}
                calculateInterpedArea(track, quantityNames{j});    
            case {'totalTime'}
                pt1 = track.pt(1);
                pt2 = track.pt(end);
                track.dq.totalTime = pt2.et - pt1.et;
            case 'mapinterpedtopts'
                mapInterpedToPts(track); %qvec(j) is the index of [track.pt.et] closets to eti(j)
            otherwise
                disp (['I don''t recognize the quantity: ' quantityNames{j}]);
        end%switch
    end%for
end %cdq
function calculateInterpedTime(track)    
    track.dq.eti = track.pt(1).et:track.dr.interpTime:track.pt(end).et;
end

function calculateInterpedLocation(track)
    pt = [track.pt];
    et = [pt.et];
    loc = [pt.loc];
    track.dq.iloc = single((interp1(et, double(loc)', track.dq.eti, 'linear'))');
end

function calculateSmoothedLocation(track)
    sigma = track.dr.smoothTime/track.dr.interpTime;
    track.dq.sloc = single(lowpass1D(track.dq.iloc, sigma));
end

function caclulateVelocity(track) 
    sigma = track.dr.derivTime/track.dr.interpTime;
    track.dq.vel = double(deriv(track.dq.sloc, sigma))/track.dr.interpTime; %velocity is in pixels per second
    track.dq.speed = double(sqrt(sum(track.dq.vel.^2, 1)));
    track.dq.nspeed = track.dq.speed / median(track.dq.speed);
    l = track.dq.speed;
    track.dq.vnorm = double(track.dq.vel./[l;l]);
    track.dq.theta = double(atan2(track.dq.vnorm(2,:), track.dq.vnorm(1,:)));
    track.dq.nvel = track.dq.vel / sqrt(mean(track.dq.speed.^2));
end

function calculateAdjSpeed(track)
    lrtime = track.dr.smoothTime * 3;
    if (~isempty(track.isrun) && any(track.isrun))
        track.dq.adjspeed = lowpass1D(track.dq.speed / median(track.dq.speed(track.isrun)), lrtime/track.dr.interpTime);
    else
        track.dq.adjspeed = lowpass1D(track.dq.nspeed, lrtime/track.dr.interpTime);
    end
end

function calculateSpeedDiff(track)
    npts = 240/track.dr.interpTime;
    if (length(track.dq.speed) < npts)
        spfilt = repmat(median(track.dq.speed), size(track.dq.speed));
    else
        spfilt = medfilt2(track.dq.speed, [1 npts], 'symmetric');
    end
    track.dq.speed_diff_local = track.dq.speed - spfilt;
    
end

    
function calculateAcceleration(track) 
    sigma = track.dr.derivTime/track.dr.interpTime;
    track.dq.deltatheta = single (deriv(unwrap(track.dq.theta), sigma))/track.dr.interpTime;
    track.dq.ddtheta = single (deriv(track.dq.deltatheta, sigma))/track.dr.interpTime;
    track.dq.acc = single(deriv(track.dq.vel,sigma))/track.dr.interpTime;
    track.dq.curv = single((track.dq.vel(1,:).*track.dq.acc(2,:) - track.dq.vel(2,:).*track.dq.acc(1,:))./(track.dq.speed.^3));
end

function calculateCovariance(track) 
    pt = [track.pt];
    et = [pt.et];
    cov = [pt.cov];
    valid = find(all(isfinite(cov), 1));
    c = double((interp1(et(valid), double(cov(:,valid))', track.dq.eti, 'linear'))');
    u = (c(1,:) + c(3,:))/2;
    d = 0.5* sqrt((c(1,:)-c(3,:)).^2 + 4 * c(2,:).^2);
    t = 0.5*asin(abs(c(2,:))./d);
    t((c(1,:) < c(3,:))) = pi/2 - t((c(1,:) < c(3,:)));
    t(c(2,:) < 0) = -t(c(2,:) < 0);
    track.dq.icov = c;
    track.dq.covMajor = sqrt(u + d);
    track.dq.covMinor = sqrt(u - d);
    track.dq.covTheta = t;
    track.dq.covRatio = track.dq.covMajor ./ track.dq.covMinor;
end
function calculateSmoothedCovariance(track) 
    sigma = track.dr.smoothTime/track.dr.interpTime;
    c = double(lowpass1D(track.dq.icov, sigma));
    u = (c(1,:) + c(3,:))/2;
    d = 0.5* sqrt((c(1,:)-c(3,:)).^2 + 4 * c(2,:).^2);
    t = 0.5*asin(abs(c(2,:))./d);
    t((c(1,:) < c(3,:))) = pi/2 - t((c(1,:) < c(3,:)));
    t(c(2,:) < 0) = -t(c(2,:) < 0);
    track.dq.scov = c;
    track.dq.scovMajor = sqrt(u + d);
    track.dq.scovMinor = sqrt(u - d);
    track.dq.scovTheta = t;
    track.dq.scovRatio = track.dq.scovMajor ./ track.dq.scovMinor;
end

function calculatePathLength(track)
   track.dq.pathLength = [0 cumsum(sqrt(sum(diff(track.dq.sloc,1,2).^2)))];
   track.dq.displacement = track.dq.sloc - repmat(track.dq.sloc(:,1), 1,length(track.dq.sloc));
end

function calculateInterpedArea(track, qn)
    if (qn(1) == 's')
        track.dq.sarea = lowpass1D(track.getDerivedQuantity('iarea'), track.dr.smoothTime/track.dr.interpTime);
    else
        pt = [track.pt];
        et = [pt.et];
        area = [pt.area];
        track.dq.iarea = single((interp1(et, double(area)', track.getDerivedQuantity('eti'), 'linear')));
    end
end

function calculatelrdtheta(track)
    lrtime = 5;
    t = track.getDerivedQuantity('theta');
    track.dq.lrdtheta = deriv(unwrap(t), lrtime/track.dr.interpTime);
end

function calculateDerivativeOfCovarianceRatio(track) 
    sc = track.getDerivedQuantity('scovRatio');
    sigma = track.dr.derivTime./track.dr.interpTime;
    track.dq.dcovRatio = deriv(sc, sigma);
end

function mapInterpedToPts(track)
    pt = [track.pt];
    et = [pt.et];
    track.calculateDerivedQuantity('eti', false);
    x = 1:length(et);
    track.dq.mapinterpedtopts = interp1(et,x,track.dq.eti, 'nearest','extrap');
end