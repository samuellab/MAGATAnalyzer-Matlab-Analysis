function calculateDerivedQuantity(track, quantityNames, recalculate)
%track.calculateDerivedQuantity(quantityName)
%calculateDerivedQuantity(track, quantityName)
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
        %short-circuit the normal track location program and use 
        %the smoothed velocity
        if strcmp(quantityNames{j}, 'vel_dp')
            calculateVel_DP(track);
            continue
        end
        if any(strcmp(quantityNames{j}, {'deltatheta', 'ddtheta', 'acc', 'curv'}))
            %this is a kludge for now to take into account that when you
            %ask for deltatheta or ddtheta, Track calculates acc & curv
            
            calculateAcceleration(track);
            %track.calculateDerivedQuantity(quantityNames{j}, false); 
            continue
        end
        if (MaggotTrack.validDQName(quantityNames{j}))
            track.calculateDerivedQuantity@MaggotTrack(quantityNames{j}, recalculate);
            continue
        end
        switch(quantityNames{j})
            case {'smoothVel', 'smoothSpeed'}
                calculateSmoothSpeed(track);
            case {'spmax', 'spmin', 'spamp'}
                calculateSpeedExtrema (track);
            otherwise
                disp (['I don''t recognize the quantity: ' quantityNames{j}]);
        end%switch
    end%for
end %cdq


function calculateSmoothSpeed(track)
    il = track.getDerivedQuantity('iloc');
    sl = lowpass1D(il, 0.5 / track.dr.interpTime);
    track.dq.smoothVel = deriv(sl, 0.25 / track.dr.interpTime)/track.dr.interpTime;
    track.dq.smoothSpeed = sqrt(sum(track.dq.smoothVel.^2));

end
function calculateVel_DP(track)
    mh = track.getDerivedQuantity('shead') - track.getDerivedQuantity('smid');
    mhnorm = sqrt(sum(mh.^2));
    try 
        mh = mh ./ [mhnorm;mhnorm];
    catch me
        disp(me.getReport);
        sum(track.dq.ihtValid) 
        size(mh)
        size(mhnorm)
        size([mhnorm;mhnorm])
        track
        track.dq
        
    end
    %track.dq.mhdir = mh;
    v = track.getDerivedQuantity('smoothVel');
    s = track.getDerivedQuantity('smoothSpeed');
    v = v ./ [s;s];
    track.dq.vel_dp = dot(v,mh);
end

    
function calculateAcceleration(track) 
    sigma = 0.25/track.dr.interpTime;
    track.calculateDerivedQuantity('smoothSpeed');
    vel = track.getDerivedQuantity('smoothVel');
    th = unwrap(atan2(vel(2,:), vel(1,:)));
    track.dq.deltatheta = single (deriv(th, sigma))/track.dr.interpTime;
    track.dq.ddtheta = single (deriv(track.dq.deltatheta, sigma))/track.dr.interpTime;
    track.dq.acc = single(deriv(track.dq.smoothVel,sigma))/track.dr.interpTime;
    track.dq.curv = single((track.dq.smoothVel(1,:).*track.dq.acc(2,:) - track.dq.smoothVel(2,:).*track.dq.acc(1,:))./(track.dq.smoothSpeed.^3));
end

function calculateSpeedExtrema(track)
    sp = track.getDerivedQuantity('speed');
    len = 1 / track.dr.interpTime;
    track.dq.spmax = ordfilt2(sp, len, ones([1 len]));
    track.dq.spmin = ordfilt2(sp, 1, ones([1 len]));
    track.dq.spamp = track.dq.spmax - track.dq.spmin;
end
