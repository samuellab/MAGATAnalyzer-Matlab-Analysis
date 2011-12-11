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
            case {'periAmp', 'periFreq', 'periTau'} %'periPhase', 'periMean'}
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
    track.dq.ihead = single((interp1(et, double(h'), track.dq.eti, 'linear','extrap')))';
    track.dq.imid = single((interp1(et, double(m'), track.dq.eti, 'linear','extrap')))';
    track.dq.itail = single((interp1(et, double(t'), track.dq.eti, 'linear','extrap')))';
   % track.dq.ihead = single((interp1(et(htv), double(h(:,htv))', track.dq.eti, 'linear')))';
   % track.dq.imid = single((interp1(et(htv), double(m(:,htv))', track.dq.eti, 'linear')))';
   % track.dq.itail = single((interp1(et(htv), double(t(:,htv))', track.dq.eti, 'linear')))';
%    track.dq.(qn) = single((interp1(et, double(loc)', track.dq.eti, 'linear'))');
    
end

function calculateSmoothedBody(track,qn)
    sigma = track.dr.smoothTime/track.dr.interpTime;
    track.dq.(qn) = single(lowpass1D(track.getDerivedQuantity(['i' qn(2:end)]), sigma));
end

function calculateVelocity(track, qn) 
    sigma = track.dr.derivTime/track.dr.interpTime;
    qn2 = qn; qn2(1) = 's';
    track.dq.(qn) = single(deriv(track.getDerivedQuantity(qn2), sigma))/track.dr.interpTime; %velocity is in pixels per second
   
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
    track.dq.iloc = single((interp1(et, double(loc)', track.dq.eti, 'linear'))');
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
    sl = squeeze(cumsum(sqrt(sum(diff(spine(:,[1 1:end], :),[],2).^2,1)),2));
    for j = 1:size(spine, 3)
        while (any(diff(sl(:,j)) == 0))
            sl(diff(sl(:,j))==0,j) = sl(diff(sl(:,j))==0,j) - 1E-8;
        end
        spine(:,:,j) = interp1(sl(:,j), spine(:,:,j)', linspace(min(sl(:,j)), max(sl(:,j)), ncp))';
    end
    spine = permute(spine,[3 1 2]);
    
    
    try
        ispine = interp1(et, spine, track.dq.eti,'linear','extrap');
        ispine = permute(ispine, [2 3 1]);
        track.dq.ispine = ispine;
    catch
        track
        size(spine)
    end
end
function calculateSpineLength(track)
    track.dq.spineLength = squeeze(sum(sqrt(sum(diff(track.dq.ispine,[],2).^2,1)),2))';
end
function calculatePeristalsis(track)

      %these peristalsis metrics are crude & will likely be replaced
      %by those developed at Janelia
      ih = track.getDerivedQuantity('ihead');
      it = track.getDerivedQuantity('itail');
      v = track.getDerivedQuantity('vnorm');
%      ht = sqrt(sum((ih-it).^2));
      th = dot(ih-it, v);
      thl = lowpass1D(th, 1);
      %derivative of the head tail distance along the velocity
      dth = deriv(th, 1);
      zcrp = find(diff(sign(dth)) > 0 & dth(1:end-1) ~= 0);
      zcrn = find(diff(sign(dth)) < 0 & dth(1:end-1) ~= 0);
      
      %interpolate to find zero crossing more precisely
      d0 = dth(zcrp);
      d1 = dth(zcrp + 1);
      zcrpi = zcrp + d0./(d0-d1);
      
      d0 = dth(zcrn);
      d1 = dth(zcrn + 1);
      zcrni = zcrn + d0./(d0-d1);
      
      ptime = interp1(track.dq.eti, zcrpi, 'linear');
      ntime = interp1(track.dq.eti, zcrni, 'linear');
      
      pper = interp1(ptime, deriv(ptime, 1), track.dq.eti, 'linear', 'extrap');
      nper = interp1(ntime, deriv(ntime, 1), track.dq.eti, 'linear', 'extrap');
      
      %period is space between zero crossings
      periTau = 0.5*(pper + nper);
      
      if (zcrp(1) < zcrn(1))
          zcrpi = zcrpi(2:end);
          zcrp = zcrp(2:end);
      end
      zcrni = zcrni(1:length(zcrpi));
      zcrn = zcrp(1:length(zcrp));
      
      %find maximum and minimum values of smoothed and unsmoothed th near zero crossing
      tm = interp1(track.dq.eti, 0.5*(zcrpi + zcrni), 'linear');
      maxth = max([thl(zcrn);thl(zcrn+1);th(zcrn);th(zcrn+1)]);
      minth = min([thl(zcrp);thl(zcrp+1);th(zcrp);th(zcrp+1)]);
      
      %amplitude is difference between maximum and minimum head tail
      %distance
      periAmp = interp1(tm, maxth-minth, track.dq.eti);
      
      sigma = median(periTau)./track.dr.interpTime;
      
      track.dq.periAmp = lowpass1D(periAmp, sigma);
      track.dq.periTau = lowpass1D(periTau, sigma);
      track.dq.periFreq = 1./track.dq.periTau;
      
%     sp = track.getDerivedQuantity('speed');
%     inds = find(abs(diff(sign(sp - median(sp)))));
%     T = 2*median(diff(track.dq.eti(inds)));
%     mysine = @(x,xdata) (x(1)*sin(2*pi*x(2)*xdata + x(3))) + x(4);
%     op = optimset('lsqcurvefit');
%     op.Display = 'off';
%     %x = lsqcurvefit(mysine, [mean(abs(sp)) 1/T 0 0], track.dq.eti, sp, [],[],op); 
%     npts = ceil(T/track.dr.interpTime);
%     x = [mean(abs(sp-median(sp))) 1/T 0 median(sp)];
%     lb = [0 0.5/T -pi 0];
%     ub = [max(sp) 1.5/T 3*pi max(sp)];
%     track.dq.periAmp = zeros(size(track.dq.eti));
%     track.dq.periFreq = zeros(size(track.dq.eti));
%     track.dq.periPhase = zeros(size(track.dq.eti));
%     track.dq.periMean = zeros(size(track.dq.eti));
%     for j = 1:length(track.dq.eti)
%         ind0 = max(j-npts, 1);
%         ind1 = min(ind0+2*npts, length(track.dq.eti));
%         ind0 = max(ind1 - 2*npts, 1);
%         inds = ind0:ind1;
%         xdata = track.dq.eti(inds) - track.dq.eti(j);
%         ydata = sp(inds);
%         x = lsqcurvefit(mysine, x, xdata, ydata, lb,ub,op);
%         x(3) = mod(x(3), 2*pi);
%         track.dq.periAmp(j) = x(1);
%         track.dq.periFreq(j) = x(2);
%         track.dq.periPhase(j) = x(3);
%         track.dq.periMean(j) = x(4);
%     end
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
     track.dq.spineTheta = track.dq.spineLength.*track.dq.spineCurv;
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