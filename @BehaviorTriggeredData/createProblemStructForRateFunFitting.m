function pd = createProblemStructForRateFunFitting (vs, timeField, polynomialDegree, trange, exprange, deltaT)
%function pd = createProblemStructForRateFunFitting (vs, timeField, polynomialDegree, trange, exprange, deltaT)
%vs = var struct as created by prepVarianceSwitchingAnalysis
ndim = length(vs);
if (ndim == 1)
    pd.ratefun = @(x, xdata) exp(polyval(x, xdata)); %[Nx1 rate function of NxD xdata given [1xK] params]
    pd.gradlogratefun =@(x,xdata) polyval(polyder(x), xdata); %NxD gradient of log of rate
    pd.hesslogratefun =@(x,xdata) polyval(polyder(polyder(x)), xdata); %NxDxD hessian of log of rate;
    pd.params_0 = zeros(1,polynomialDegree+1);
    pd.params_0(end) = 1;
else
    switch (polynomialDegree)
        case 1
            pd.ratefun = @(x, xdata) exp(xdata * x(1:(end-1))' + x(end)); %[Nx1 rate function of NxD xdata given [1xK] params]
            pd.gradlogratefun = @(x,xdata) repmat(x(1:(end-1)), [size(xdata,1), 1]); %NxD gradient of log of rate
            pd.hesslogratefun = @(x,xdata) zeros([size(xdata,1) size(xdata,2) size(xdata,2)]); %NxDxD hessian of log of rate;
            pd.params_0 = zeros(1,ndim+1);
            for j = 1:ndim
                pd.params_0(j) = mean(vs(j).turn.x_conv)/var(vs(j).noturn.x_conv);
            end
            pd.params_0(end) = 1;
        case 2
            if (ndim == 2)
                pd.ratefun = @(x,xdata) exp(x(1)*xdata(:,1).^2 + x(2)*xdata(:,1).*xdata(:,2) + x(3)*xdata(:,2).^2 + x(4)*xdata(:,1) + x(5)*xdata(:,2) + x(6));
                pd.gradlogratefun = @(x,xdata) [2*x(1)*xdata(:,1) + x(2)*xdata(:,2) + x(4);x(2)*xdata(:,2) + 2*x(3)*xdata(:,2) + x(5)];
                pd.hesslogratefun = @twobytwohess; 
                pd.params_0 = zeros(1,6);
                pd.params_0(end) = 1;
            else
                error ('polynomial of degree > 1 not supported for data of dimenion > 2');
            end
        otherwise
            error ('polynomial of degree > 2 not supported for data of dimension > 1');
    end
end
pd.timeField = timeField;
pd.temporalratemod = @(x,tdata) x(1)*exp(-tdata/x(2)) + 1;

turn = [vs.turn];
run = [vs.noturn];
if (length(turn) > 1) && any(any(diff([turn.eti],[],2))) || any(any(diff([run.eti], [], 2)))
    error ('I expect that data in the different var structures should be temporally aligned');
end

pd.turnT = turn(1).(timeField);
pd.turnEti = turn(1).eti;
pd.turnX = [turn.x_conv];
pd.turnExpnum = turn(1).expnum;

pd.runT = run(1).(timeField);
pd.runEti = run(1).eti;
pd.runX = [run.x_conv];
pd.runExpnum = run(1).expnum;

existsAndDefault('deltaT', median(diff(pd.runEti)));
pd.deltaT = deltaT;
%pd.deltaT = median(diff(pd.runEti));

period = vs(1).period;

if (strcmpi(timeField, 'ton') || strcmpi(timeField, 'toff'))
   pd.tx = 0:pd.deltaT:period;   
   pd.pad = true;
   if (~exist('trange', 'var') || isempty('trange'))
       minTime = min(run(1).eti(run(1).ton >= 0 & run(1).toff >= 0));
       trange = minTime:period:max(run(1).eti + pd.deltaT); trange = trange([1 end]);
   end
    trange(1) = max(trange(1),min(pd.runEti));
    trange(2) = min(trange(2),max(pd.runEti));
else
    existsAndDefault('trange', [max(vs(1).period - vs(1).tshift,min(pd.runEti)) max(pd.runEti)+pd.deltaT]);
    trange(1) = max(trange(1),min(pd.runEti));
    trange(2) = min(trange(2),max(pd.runEti));
    pd.tx = (trange(1):pd.deltaT:trange(end))';
    pd.pad = false;
end
existsAndDefault ('exprange', [min(pd.runExpnum) max(pd.runExpnum)]);
pd.exprange = exprange;

pd.timeRange = trange;
pd.etx = linspace(min(trange), max(trange), 100);

ti = pd.turnEti >= min(trange) & pd.turnEti < max(trange) & pd.turnExpnum >= min(exprange) & pd.turnExpnum <= max(exprange);
pd.turnT = pd.turnT(ti);
pd.turnEti = pd.turnEti(ti);
pd.turnX = pd.turnX(ti,:);
pd.turnExpnum = pd.turnExpnum(ti);

ri = pd.runEti >= min(trange) & pd.runEti < max(trange) & pd.runExpnum >= min(exprange) & pd.runExpnum <= max(exprange);
pd.runT = pd.runT(ri);
pd.runEti = pd.runEti(ri);
pd.runX = pd.runX(ri,:);
pd.runExpnum = pd.runExpnum(ri);

%pd.runX = pd.runX - mean(pd.runX);
%pd.turnX = pd.turnX - mean(pd.turnX);

pd.params_0(end) = log(length(pd.turnT)/(length(pd.runT)*pd.deltaT));

pd.Q_alpha = 0.2*eye(ndim);
pd.Q_s = pd.Q_alpha;

pd.tparams_0 = [0 1000];
pd.alpha_0 = ones(length(pd.tx), ndim);
pd.sigma_0 = 10*ones(length(pd.tx), ndim);
pd.v0 = 1;
pd.w_0 = pd.Q_alpha*1000*eye(ndim);
pd.maxreps = 3;

pd.separateExperiments = false;
pd.period = vs(1).period;
pd.tshift = vs(1).tshift;

function hess = twobytwohess (x, xdata)
hess = zeros(size(xdata,1), 2, 2);
hess(:,1,1) = 2*x(1);
hess(:,1,2) = x(2);
hess(:,2,1) = x(2);
hess(:,2,2) = x(3);