function calculateDerivedQuantity(track, quantityNames, recalculate)
% calculate a quantity based on extracted position and posture
% track.calculateDerivedQuantity(quantityName, recalculate)
% calculateDerivedQuantity(track, quantityName, recalculate)
%
% calculates a derived quantity (for a list of derived quantities, use
% MaggotTrack/validDQName)
% if the quantity has already been calculated (track.dq.(quantityName)
% exists), then the value is not recalculated, unless recalculate is true
% "recalculate" may or may not cause other quantities to also be
% recalculated
%
% outputs: none
% inputs: 
%       TRACK < MaggotTrack
%       QUANTITYNAMES: string naming quantity, or cell of strings naming
%           quantitites
%       RECALCULATE < bool


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
        %the midpoint instead of the center of mass
        if strcmp(quantityNames{j}, 'iloc')
            calculatePosition(track);
            continue
        end
        if (Track.validDQName(quantityNames{j}))
            track.calculateDerivedQuantity@Track(quantityNames{j}, recalculate);
            continue
        end
        switch(quantityNames{j})
            case {'ihead', 'itail', 'imid'}
                track.calculateDerivedQuantity({'eti', 'ihtValid'}, false);
                calculateInterpedBody(track);
            case {'shead', 'stail', 'smid'}
                qname = quantityNames{j};
                qname(1) = 'i';
                track.calculateDerivedQuantity(qname, false);
                calculateSmoothedBody(track, quantityNames{j});
            case {'vhead', 'vtail', 'vmid'}
                qname = quantityNames{j};
                qname(1) = 'i';
                track.calculateDerivedQuantity(qname, false);
                calculateVelocity(track, quantityNames{j});
            case {'sphead', 'sptail', 'spmid'}
                calculateSpeed(track, quantityNames{j});
            case {'ibodytheta','sbodytheta'}
                track.calculateDerivedQuantity({'ihead', 'itail', 'imid'});
                if (quantityNames{j}(1)=='s' || quantityNames{j}(1) == 't')
                    track.calculateDerivedQuantity({'shead', 'stail', 'smid'});
                end
                calculateBodyAngle(track, quantityNames{j})
            case {'dbodytheta', 'dsbodytheta'}
                calculateDBodyAngle(track);
            case {'ihtValid'}
                track.calculateDerivedQuantity('eti');
                calculateInterpedHTValid(track);
            case {'vel_dp'}
                track.calculateDerivedQuantity({'shead', 'smid', 'vel', 'speed'});
                calculateVel_DP(track);
            case {'vheadperp','spheadperp'}
                track.calculateDerivedQuantity({'shead', 'smid', 'stail'});
                calculateVHeadPerp(track);
            case {'itmdir', 'stmdir', 'imhdir', 'smhdir'}
                track.calculateDerivedQuantity({'ihead', 'itail', 'imid'});
                if (quantityNames{j}(1)=='s' || quantityNames{j}(1) == 't')
                    track.calculateDerivedQuantity({'shead', 'stail', 'smid'});
                end
                calculateBodyDirection(track, quantityNames{j});
            case 'ispine'
                track.calculateDerivedQuantity('eti');
                calculateInterpedSpine(track);
            case 'spineLength'
                track.calculateDerivedQuantity({'eti', 'ispine'});
                calculateSpineLength(track);
            case {'spineCurv'}
                track.calculateDerivedQuantity({'eti', 'ispine', 'spineLength'});
                calculateSpineCurv(track);
            case {'spineTheta'}
                calculateSpineTheta(track);
            case {'sspineTheta'}
                track.dq.sspineTheta = lowpass1D(track.getDerivedQuantity('spineTheta'), track.dr.smoothTime/track.dr.interpTime);
                
            case {'dspineCurv', 'dspineTheta'}
                calculateDSpineCurv(track);
                
            case 'spineWidth'
                track.dq.spineWidth = track.getDerivedQuantity('iarea')./track.getDerivedQuantity('spineLength');
            case {'fastTailSpeed', 'fastHeadSpeed', 'fastMidSpeed'}
                calculateFastSpeed(track);
            case {'periAmp', 'periFreq', 'periPhase'} %'periPhase', 'periTau', 'periMean'}
                calculatePeristalsis(track);
            case {'spineDist'}
                calculateSpineDist (track);
            
            otherwise
                disp (['I don''t recognize the quantity: ' quantityNames{j}]);
        end%switch
    end%for
end %cdq
function calculateInterpedHTValid(track)    
    pt = [track.pt];
    et = [pt.et];
    htv = double([pt.htValid]);
    
    %et = [track.pt.et];
    %htv = double([track.pt.htValid]);
    track.dq.ihtValid =  ((interp1(et, htv', track.dq.eti, 'linear')) > 0.99);
end

function calculateInterpedBody(track)
    pt = [track.pt];
    htv = [pt.htValid];
    pt = pt(htv);
    et = [pt.et];
    h = [pt.head];
    t = [pt.tail];
    m = [pt.mid];
    %et = [track.pt([track.pt.htValid]).et];
    %loc = [track.pt([track.pt.htValid]).(qn(2:end))];
    if (length(et) < 2)
        track.dq.ihead = NaN(2,length(track.dq.eti));
        track.dq.imid = NaN(2,length(track.dq.eti));
        track.dq.itail = NaN(2,length(track.dq.eti));
        return;
    end
    try
        track.dq.ihead = double((interp1(et, double(h'), track.dq.eti, 'linear','extrap')))';
        track.dq.imid = double((interp1(et, double(m'), track.dq.eti, 'linear','extrap')))';
        track.dq.itail = double((interp1(et, double(t'), track.dq.eti, 'linear','extrap')))';
    catch me
        disp (me.getReport());
        size(et)
        size(h)
        track.dq.ihead = NaN(2,length(track.dq.eti));
        track.dq.imid = NaN(2,length(track.dq.eti));
        track.dq.itail = NaN(2,length(track.dq.eti));
    end
   % track.dq.ihead = single((interp1(et(htv), double(h(:,htv))', track.dq.eti, 'linear')))';
   % track.dq.imid = single((interp1(et(htv), double(m(:,htv))', track.dq.eti, 'linear')))';
   % track.dq.itail = single((interp1(et(htv), double(t(:,htv))', track.dq.eti, 'linear')))';
%    track.dq.(qn) = single((interp1(et, double(loc)', track.dq.eti, 'linear'))');
    
end

function calculateSmoothedBody(track,qn)
    sigma = track.dr.smoothTime/track.dr.interpTime;
    track.dq.(qn) = double(lowpass1D(track.getDerivedQuantity(['i' qn(2:end)]), sigma));
end

function calculateVelocity(track, qn) 
    sigma = track.dr.derivTime/track.dr.interpTime;
    qn2 = qn; qn2(1) = 's';
    track.dq.(qn) = double(deriv(track.getDerivedQuantity(qn2), sigma))/track.dr.interpTime; %velocity is in pixels per second
   
end
function calculateSpeed(track, qn) 
    qn2 = qn(2:end); qn2(1) = 'v';
    track.dq.(qn) = sqrt (sum(track.getDerivedQuantity(qn2).^2));
end

    
function calculatePosition(track)
    track.calculateDerivedQuantity('eti');
    pt = [track.pt];
    et = [pt.et];
    loc = [pt.mid];
    track.dq.iloc = double((interp1(et, double(loc)', track.dq.eti, 'linear'))');
end

    
function calculateBodyAngle(track, qn)
    mh = track.dq.([qn(1) 'head']) - track.dq.([qn(1) 'mid']);
    tm = track.dq.([qn(1) 'mid'])- track.dq.([qn(1) 'tail']); 
    alltheta = [atan2(mh(2,:), mh(1,:));atan2(tm(2,:), tm(1,:))]; %tail angle - head angle
    track.dq.(qn) = -diff(unwrap(alltheta));
end

 function calculateBodyDirection(track, qn)
    mh = track.dq.([qn(1) 'head']) - track.dq.([qn(1) 'mid']);
    tm = track.dq.([qn(1) 'mid'])- track.dq.([qn(1) 'tail']); 
    if (qn(2) == 't')
        l = sqrt(sum(tm.^2));
        track.dq.(qn) = tm./[l;l];
    else
        l = sqrt(sum(mh.^2));
        track.dq.(qn) = mh./[l;l];
    end
 end


function calculateVel_DP(track)
    mh = track.dq.shead - track.dq.smid;
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
    v = track.dq.vel;
    s = track.dq.speed;
    v = v ./ [s;s];
    track.dq.vel_dp = dot(v,mh);
end
function calculateVHeadPerp(track) 
%{
    tm = track.dq.smid - track.dq.stail;
    tmnorm = sqrt(sum(tm.^2));
    vh = deriv(track.dq.shead, track.dr.derivTime./track.dr.interpTime);
    vhnorm = sqrt(sum(vh.^2));
    vhnorm(vhnorm <= eps) = eps;
    try 
        tm = tm ./ [tmnorm;tmnorm];
        vhn = vh ./ [vhnorm;vhnorm];
    catch me
        disp(me.getReport);
        sum(track.dq.ihtValid) 
        size(tm)
        size(tmnorm)
        size([tmnorm;tmnorm])
        track
        track.dq
        
    end
    alpha = dot(vhn,tm);
%}
    tm = track.getDerivedQuantity('stmdir');
    vh = track.getDerivedQuantity('vhead');
    vhnorm = sqrt(sum(vh.^2));
    vhnorm(vhnorm <= eps) = eps;
    try         
        vhn = vh ./ [vhnorm;vhnorm];
    catch me
        disp(me.getReport);
        sum(track.dq.ihtValid) 
        size(vh)
        size(vhnorm)
        size([vhnorm;vhnorm])
        track
        track.dq       
    end
    alpha = dot(vhn,tm);
    track.dq.vheadperp = vh - vh.*[alpha;alpha];
    track.dq.spheadperp = sqrt(sum(vh.^2));
end
function calculateInterpedSpine(track)    
    
    ncp = 11;
    pt = [track.pt];
    et = [pt.et];
    htv = [pt.htValid];
    pt = pt(htv);
    et = et(htv);
    
    spine = reshape([pt.spine], 2, ncp, []);
    valid = squeeze(all(all(isfinite(spine))));
    spine = spine(:,:,valid);
    et = et(valid);
    
    sl = squeeze(cumsum(sqrt(sum(diff(spine(:,[1 1:end], :),[],2).^2,1)),2));
    for j = 1:size(spine, 3)
        while (any(diff(sl(:,j)) == 0))
            sl(diff(sl(:,j))==0,j) = sl(diff(sl(:,j))==0,j) - 1E-8;
        end
        spine(:,:,j) = interp1(sl(:,j), spine(:,:,j)', linspace(min(sl(:,j)), max(sl(:,j)), ncp))';
    end
    spine = permute(spine,[3 1 2]);
    
    
   % try
        ispine = interp1(et, spine, track.dq.eti,'linear','extrap'); %use normal linear interpolation for valid spine values
        %use nearest neighbor interpolation for invalid spine values --
        %added by mhg 7/2/2012
        ihtinv = ~track.getDerivedQuantity('ihtValid');
        if (any(ihtinv))
            ispine(ihtinv, :, :) = interp1(et, spine, track.dq.eti(ihtinv), 'nearest', 'extrap');
        end    
        ispine = permute(ispine, [2 3 1]);
        
        
        
        
        track.dq.ispine = ispine;
%     catch me
%         disp(me.getReport)
%         track
%         size(spine)
%     end
end
function calculateSpineLength(track)
    track.dq.spineLength = squeeze(sum(sqrt(sum(diff(track.dq.ispine,[],2).^2,1)),2))';
end
function calculateFastSpeed(track)
  track.calculateDerivedQuantity({'eti', 'itail', 'ihead', 'imid'});
  tm = track.dq.imid - track.dq.itail;
  tm = tm./repmat(sqrt(sum(tm.^2)),[2 1]);
  %mh = track.dq.ihead - track.dq.imid;
  %th = track.dq.ihead - track.dq.itail;
  track.dq.fastTailSpeed = dot (deriv(track.dq.itail, 1), tm)/track.dr.interpTime;
  track.dq.fastHeadSpeed = dot (deriv(track.dq.ihead, 1), tm)/track.dr.interpTime;
  track.dq.fastMidSpeed = dot (deriv(track.dq.imid, 1), tm)/track.dr.interpTime;
  
  %sqrt(sum(deriv(track.dq.itail,1)).^2)/track.dr.interpTime;
  %track.dq.fastHeadSpeed = sqrt(sum(deriv(track.dq.ihead,1)).^2)/track.dr.interpTime;
  %track.dq.fastMidSpeed = sqrt(sum(deriv(track.dq.imid,1)).^2)/track.dr.interpTime;
end
function calculatePeristalsis(track)

      %tail velocity seems smoothest
      track.calculateDerivedQuantity({'eti', 'itail'});
      derivtime = 0.05;
      vt = (deriv(track.dq.itail,derivtime/track.dr.interpTime))/track.dr.interpTime;
      tm = lowpass1D(track.dq.imid - track.dq.itail, derivtime/track.dr.interpTime);
      
      vtm = dot(vt, tm);
%      vtm = vt;% - mean(vt);
      vtm(vtm > percentile(vtm, 0.99)) = percentile(vtm, 0.99);
      vtm(vtm < percentile(vtm, 0.01)) = percentile(vtm, 0.01);
      time_window = min(track.dq.eti(end)-track.dq.eti(1), 10);
      
      %using all tail velocity, find overall best peristaltic frequency
      Hs = spectrum.welch('Hamming', time_window/track.dr.interpTime);
      hpsd = Hs.psd(vtm, 'Fs', 1 / track.dr.interpTime, 'NormalizedFrequency', false);
      ps = hpsd.Data;
      f = hpsd.Frequencies;
      ps = ps(f > 0.5);
      f = f(f > 0.5);
      [~,I] = max(ps); 
      bestfreq = f(I);
      if (bestfreq < 1 || bestfreq > 5)
          warning (['overall peristaltic frequency: ' num2str(bestfreq,3) ' is outside expected range']);
      end
      %now use autocorrelation to find local frequency
      periodInPts = 1/(track.dr.interpTime*bestfreq);
      minshift = floor(periodInPts * 0.8);
      maxshift = ceil(periodInPts * 1.4);
      winsize = periodInPts*1.5;
      acf = autocorrInWindow(vtm, minshift, maxshift, winsize, true);
      ssf = squareSumInWindow(vtm, minshift, maxshift, winsize, true);
      z = 2*acf./ssf;
      [b,I] = max(z);
      a = interp2(z, 1:length(I), I-1, '*nearest');
      c = interp2(z, 1:length(I), I+1, '*nearest');
      a(isnan(a)) = b(isnan(a));
      c(isnan(c)) = b(isnan(c));
      I2 = I + 0.5 * (a-c)./(a-2*b+c);
      fs = 1./(track.dr.interpTime*(minshift - 1 + I2));
      zcut = 0.75;%percentile(b, 0.1);
%       figure(1)
%       pcolor (track.dq.eti, minshift:maxshift, z); shading flat; colorbar vert;
%       figure(2)
%       pcolor (track.dq.eti, minshift:maxshift, sort(z)); shading flat; colorbar vert;
%       zz = sort(z, 1 , 'descend');
%       figure(3); clf
%       for n = 2:10
%           plot (track.dq.eti, mean(zz(1:n,:))); hold all
%       end
      pv = b > zcut & (max(z) > 2*min(z));
      if (nnz(pv) < 2)
          pv = b > zcut; %this is a kludge -- not sure what pv, b, zcut are 10/26/2014 -- gershow
      end
      
      pf = interp1(track.dq.eti(pv) + (winsize/2.0)*track.dr.interpTime, fs(pv), track.dq.eti, 'linear', 'extrap');
      targetFreq = 1/50;
      
      %stretch time to create a uniform frequency
      dfun = @(xx,ss) targetFreq/track.dr.interpTime./interp1(pf, ss, 'linear', bestfreq);
      nspts = bestfreq*track.dr.interpTime/targetFreq * length(pf);
      [~,s] = ode15s(dfun, 1:(nspts*1.25), 1);
      s = s(s < length(pf))';
      
      %tail velocity and position semi-regularized to have a peristaltic frequency of
      %targetFreq regardless of position
      vti = interp1(vtm, s, 'linear', 0);
      %tpi = interp1(track.dq.itail', s, 'linear')';
      
      %find period of peristaltic movement when tail is relatively still
      tk = unique([0:(0.5/targetFreq) -(0:(0.5/targetFreq))]);
      ck = exp(sqrt(-1)*2*pi*targetFreq*tk);
      c = conv2(vti, ck, 'same');
      th = angle(c);
      %mag = abs(c);
      [~,fi] = findPositiveZeroCrossings(th);
      periInds = interp1(s, fi); %points in non-morphed time that mark pause point of peristalsis cycle
      periDist = sqrt(sum(diff(interp1(track.dq.imid', periInds))'.^2));
      periPeriod = diff(interp1(track.dq.eti, periInds));
      periQuality = interp1(b, periInds);
      periQuality = interp1(periQuality, 1.5:length(periQuality), 'linear');
      periMag = interp1(interp1(abs(c), fi),1.5:length(periInds)); 
      
      
      periTimePoint = 0.5 * (interp1(track.dq.eti, periInds(1:(end-1))) + interp1(track.dq.eti, periInds(2:end)));
      track.dq.periAmp = interp1(periTimePoint, periDist, track.dq.eti, 'linear', 'extrap');
      track.dq.periPeriod = interp1(periTimePoint, periPeriod, track.dq.eti, 'linear', 'extrap');
      track.dq.periFreq = pf;
      track.dq.periValid = pv;
      %figure(1); clf; plot (pf, track.dq.periFreq, 'b.');
      
      track.dq.periPhase = interp1(s, th, 1:length(track.dq.eti));
      track.dq.periSpeed = interp1(periTimePoint, periDist./periPeriod, track.dq.eti, 'linear', 'extrap');
      track.dq.periQuality = b;
%       figure(1);
%       plot (periPeriod, periDist, 'b.');
%       figure(2);
%       plot (periPeriod, periQuality, 'r.');
%       figure(3);
%       plot (periQuality, periDist, 'g.');
%       figure(4)
%       plot (periPeriod(periQuality > zcut), periDist(periQuality > zcut), 'b.');
       
      
      
      
end

function calculateSpineDist (track)
    is = track.getDerivedQuantity('ispine');
    is = is - repmat(permute(track.getDerivedQuantity('iloc'), [1 3 2]), [1 size(is,2) 1]);
    splen = size(is, 2);
    ds = diff(is, [],3);
%    dmp = diff(track.getDerivedQuantity('iloc'),[],2);
    
    dsq = sqrt(sum(ds.^2));
 %   mpdist = sqrt(sum(diff(track.getDerivedQuantity('iloc'),[],2).^2)); 
    sd = (squeeze(sum(dsq,2)))';
  %  track.dq.spineDist = ([0 (sd/splen - mpdist)])/track.dr.interpTime;
  track.dq.spineDist = ([0 sd/(splen)]);
end

function calculateDBodyAngle(track) 
    sigma = track.dr.derivTime/track.dr.interpTime;
    th = unwrap(track.getDerivedQuantity('ibodytheta'));
    track.dq.dbodytheta = deriv(th, sigma);
    th = unwrap(track.getDerivedQuantity('sbodytheta'));
    track.dq.dsbodytheta = deriv(th, sigma);
    

end

function calculateSpineCurv(track)
     is = track.getDerivedQuantity('ispine');
     sz = size(is,2);
     is = interp1(permute(is, [2 1 3]), linspace(1,sz,100));
     v = 0.5*(diff(is(2:end,:,:)) + diff(is(1:end-1,:,:)));
     a = 0.5*(diff(v(2:end,:,:)) + diff(v(1:end-1,:,:)));
     v = v(2:end-1,:,:);
     cv = mean((v(:,1,:).*a(:,2,:)) - (v(:,2,:).*a(:,1,:))./(sum(v.^2, 2)).^1.5);
     track.dq.spineCurv = squeeze(cv)';
 %    track.dq.spineTheta = track.dq.spineLength.*track.dq.spineCurv;
end

function calculateSpineTheta(track)
    is = track.getDerivedQuantity('ispine');
    nspinepts = size(is,2);

    range = max(2,round(nspinepts/5)); 
    midind = ceil(nspinepts/2) + (-range:range);

    sqe = zeros(length(midind), size(is, 3));
    dt = sqe;
    %thh = sqe;
    for j = 1:length(midind)

        tail = is(:,1:midind(j), :);
        head = is(:,midind(j):end, :);
        [~, dvt, sqet] = fitLine (tail);
        [~, dvh, sqeh] = fitLine (head);
        tht = atan2(dvt(2,:), dvt(1,:));
        thh = atan2(dvh(2,:), dvh(1,:));
        dt(j,:) = diff(unwrap([tht;thh])); 
        sqe(j,:) = sqet+sqeh;
    end

    [~,I] = min(sqe);
    track.dq.spineTheta = dt(sub2ind(size(dt), I, 1:length(I)));
end

function calculateDSpineCurv(track)
    sigma = track.dr.derivTime/track.dr.interpTime;
    th = unwrap(track.getDerivedQuantity('spineTheta'));
    track.dq.spineTheta = deriv(th, sigma);
    cv = track.getDerivedQuantity('spineCurv');
    track.dq.dspineCurv = deriv(cv, sigma);
    
end