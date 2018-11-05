

clear zz zz_big zz_42a ZZ zz_berlin Z_42a Z_berlin dts taus dts2 taus2
new_axes = 0; % if you want to use other dt,tau axes and interpolate

%LLs = [LogLBerlin, LogL42a_fine];
LLs = [LogL_Berlin, LogL_42a];

for k=1:length(LLs)

LLstruct = LLs(k);
dts{k} = LLstruct.dts;
taus{k} = LLstruct.taus;

nll{k} = LLstruct.LL.logP_fitR;
zz{k} =  -nll{k} + min(min(nll{k}));

a1 = [.05 .1:.1:.3];
a3 = [1.2 2 3];
a2 = .5:.1:1;
dts_new = [a2];
taus_new = taus{k};

[lxx, lxy] = meshgrid(taus{k}, dts{k});
[lxx2, lxy2] = meshgrid(taus_new, dts_new);
lxdata = [lxx(:) lxy(:)];

lxdata2 = [lxx2(:) lxy2(:)];
nll_new{k} = interp2(lxx, lxy, nll{k}, lxx2, lxy2, 'spline');
zz_new{k} =  -nll_new{k} + min(min(nll_new{k}));
end

if(new_axes)
    zz_42a = zz_new{2};
    zz_berlin = zz_new{1};
    tauaxis_42a = taus_new;
    tauaxis_ber = taus_new;
    dtaxis_42a = dts_new;
    dtaxis_ber = dts_new;
    ldata = lxdata2;
    ldata2 = lxx2;
else
    zz_42a = zz{2};
    zz_berlin = zz{1};
    tauaxis_42a = taus{2};
    dtaxis_42a = dts{2};
    tauaxis_ber = taus{1};
    dtaxis_ber = dts{1};
    ldata = lxdata;
    ldata2 = lxx;
end
    

Z_42a = zz_42a;
Z_berlin = zz_berlin;

% fit log-likelihood to a quadratic (to the log of a 2-dim normal distribution)
% fitfun == C + log( normPDF(x,y) )

fitfun = @(x, xd) x(6) - log(2*pi*sqrt(x(3)^2*x(4)^2 - x(5)^2)) - (x(3)^2*x(4)^2)/(2*x(3)^2*x(4)^2-x(5)^2) * ( (xd(:,1)-x(1)).^2./x(3)^3 -2*x(5)*(xd(:,1)-x(1)).*(xd(:,2)-x(2))./(x(3)^2*x(4)^2) + (xd(:,2)-x(2)).^2./x(4)^2 );
op = optimoptions('lsqcurvefit');
op.MaxFunctionEvaluations = 1e4;
op.MaxIter = 1e4;

problem2.solver = 'lsqcurvefit';
problem2.options = op;
problem2.objective = fitfun;
problem2.ub = [20 2 20 5 10 Inf];
problem2.lb = [0 0 0 0 0 -Inf];

% berlin: trim logL to only fit the region where logL>=-7 (relative to best value)
Z_berlin_fine = Z_berlin(6:8,4:9);
dts_ber = dtaxis_ber(6:8);
taus_ber = tauaxis_ber(4:9);

[lxx, lxy] = meshgrid(taus_ber, dts_ber);
lxdata = [lxx(:) lxy(:)];

x0 = [8 .5 5 .1 0.1 min(min(Z_berlin_fine))];
problem2.xdata = lxdata;
problem2.ydata = Z_berlin_fine(:);
problem2.x0 = x0;
fitP_berlin = lsqcurvefit(problem2);

% 42a
Z_42a_fine = Z_42a(:,:);
dts_42a = dtaxis_42a(:);
taus_42a = tauaxis_42a(:);


[lxx, lxy] = meshgrid(taus_42a, dts_42a);
lxdata = [lxx(:) lxy(:)];

x0 = [8 .5 5 .1 0.1 min(min(Z_42a_fine))];
problem2.xdata = lxdata;
problem2.ydata = Z_42a_fine(:);
problem2.x0 = x0;
fitP_42a = lsqcurvefit(problem2);



% covariance matrix

cm_ber = [round(fitP_berlin(3))^2 fitP_berlin(5); fitP_berlin(5) fitP_berlin(4)^2]; % covariance matrix
[V, D] = eig(cm_ber);
D_ber = abs(D);

cm_42a = [round(fitP_42a(3))^2 fitP_42a(5); fitP_42a(5) fitP_42a(4)^2]; % covariance matrix
[V, D] = eig(cm_42a);
D_42a = abs(D);

ss = 4.605; % 90% confidence
ss = 5.991; % 95% confidence

% confidence regions: see http://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
conf_berlin_dt = 2*sqrt(ss*D_ber(1,1));
conf_berlin_tau = 2*sqrt(ss*D_ber(2,2));

conf_42a_dt = 2*sqrt(ss*D_42a(1,1));
conf_42a_tau = 2*sqrt(ss*D_42a(2,2));


return
