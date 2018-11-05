function calculateMetrics(st)
%WormSharpTurn/calculateMetrics(st)
%
%calculates the properties of the sharp turn;
track = st.track;
track.calculateDerivedQuantity({'scovRatio','sloc','theta','deltatheta'});
wso = track.so;

 
s = max(st.startInd-wso.ptBuffer, 1);
e = min(st.endInd+wso.ptBuffer, length(track.dq.sloc));


st.thetaIn = track.dq.theta(s);
st.thetaOut = track.dq.theta(e);
st.dTheta = diff(unwrap([st.thetaIn;st.thetaOut]));
  

[~,lm] = max(abs(track.dq.deltatheta(st.startInd:st.endInd)));
st.centralInd = st.startInd - 1 + lm;
st.inds = st.startInd:st.endInd;


%old: flag on total angle change
st.typeCode = -1;
if (abs(st.dTheta) < wso.alignedTheta && (min(track.dq.scovRatio(st.startInd:st.endInd)) > wso.reversalCovThresh) )
    st.typeCode = 0;
end

if (abs(st.dTheta) > pi - wso.alignedTheta && min(track.dq.scovRatio(st.inds)) > wso.reversalCovThresh)
    st.typeCode = 1;
end

if (isfield(wso, 'omegaCovThresh') && ~isempty(wso.omegaCovThresh))
    ct = wso.omegaCovThresh;
else
    ct = wso.reversalCovThresh;
end

if (min(track.dq.scovRatio(st.inds)) > ct)
    if (abs(st.dTheta) < pi / 2)
        st.typeCode = 0;
    else
        st.typeCode = 1;
    end
end

if (st.typeCode == -1)
    touchUpOmegaTurn(st);
end
    
if (st.typeCode == 1)
    touchUpReversal(st);
end

st.loc = track.dq.sloc(:,st.centralInd);
st.inds = st.startInd:st.endInd;
st.dTheta = diff(unwrap([st.thetaIn;st.thetaOut]));



function touchUpOmegaTurn (st)
    
cr = st.track.getDerivedQuantity('scovRatio');
dcr = st.track.getDerivedQuantity('dcovRatio');
[mcr,I] = min(cr(st.inds));
if (isempty(I))
    st
    st.track
    pause
end
%cthresh = 0.5 * mcr + 0.25*(cr(st.startInd) + cr(st.endInd));
cI = I + st.startInd - 1;

st.centralInd = cI;
start = find(dcr(1:cI) > 0 & cr(1:cI) > st.track.so.reversalCovThresh & cr(1:cI) > cr(cI), 1, 'last');
if (isempty(start))
    start = st.startInd;
end
stop = find(dcr(cI:end) < 0 & cr(cI:end) > st.track.so.reversalCovThresh & cr(cI:end) > cr(cI), 1, 'first');

if (isempty(stop))
    stop = st.endInd;
else
    stop = cI + stop - 1;
end

st.startInd = start;
st.endInd = stop;
st.inds = start:stop;
    
st.thetaIn = st.track.dq.theta(start);
st.thetaOut = st.track.dq.theta(stop);

function touchUpReversal(st)
    
dt = st.track.getDerivedQuantity('deltatheta');
s = sign(dt(st.centralInd));
ddt = s*st.track.getDerivedQuantity('ddtheta');

si = find(sign(dt(1:(st.centralInd-1))) ~= s | ddt(1:(st.centralInd-1)) < 0, 1, 'last');
if (~isempty(si))
    st.startInd = si + 1;
end

ei = find(sign(dt((st.centralInd + 1):end)) ~= s | ddt((st.centralInd + 1):end) > 0, 1, 'first');
if (~isempty(ei))
    st.endInd = st.centralInd + ei;
end
