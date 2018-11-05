function [theta, w, tx, diagnostics] = estimateTimeVaryingROG (turnTime, turnX, runTime, runX, q_rate, samplingInterval, tx, initTheta, initW)
%function [theta, w] = estimateTimeVaryingROG (turnTime, turnX, runTime, runX, q_rate, samplingInterval, tx)
%
%implements adaptive point process filter from Eden, U.T., Frank, L.M., Barbieri, R., Solo, V., and Brown, E.N. (2004). Dynamic Analysis of Neural Encoding by Point Process Adaptive Filtering. Neural Computation 16, 971–998.
%for ratio of gaussians function:
%r1(x) = exp(theta(3) + 0.5 * x^2 - 0.5*theta(1)*(x-theta(2))^2 (rateCalcs)
%or r2(x) = exp(theta(3) + theta(2)*x + theta(1)*x^2) (rate2Calcs)
%or r3(x1, x2) = exp( theta(1)*x1^2 + theta(2)*x2^2 + theta(3)*x1*x2 +%theta(4)*x2 + theta(5)*x2 + theta(6) )  if size(turnX)>2

if size(turnX, 2) > 1 %for two (2) white-noise stimuli
    
    turnTime = turnTime';
    turnX = turnX';
    runTime = runTime';
    runX = runX';
    allTime = [turnTime runTime];
    allX(1,:) = [turnX(1,:) runX(1,:)];
    allX(2,:) = [turnX(2,:) runX(2,:)];
    
    if (nargin < 6)
        samplingInterval = median(diff(allTime));
    end
    if (nargin < 7)
        tx =min(allTime):samplingInterval:max(allTime);
    end
    dt = diff(binEdgesFromCenters(tx));
    
    [~,~,turnBin] = histcounts(turnTime,  binEdgesFromCenters(tx));
    [~,~,allBin] = histcounts(allTime,  binEdgesFromCenters(tx));
    
    w = zeros(6,6,length(dt));
    if (~exist('initW', 'var'))
        w(:,:,1) = diag([1 1 1 1 1 1]);
    else
        w(:,:,1) = initW;
    end
    theta = zeros(6,length(dt));
    if (~exist('initTheta', 'var'))
        r0 = log(length(turnTime)) - log(length(allTime)*samplingInterval);
        theta(:,1) = [0;0;0;0;0;r0];
    else
        theta(:,1) = initTheta;
    end
    
    for j = 1:length(dt)
        k = max(j-1,1);
        turninds = turnBin == j;
        allinds = allBin == j;
        %theta_stim(j) = atan2(mean(turnX(2, turninds)), mean(turnX(1, turninds)));
        [lambdaT, dlambdaT] = Combined_rateCalcs(turnX(1,turninds), turnX(2,turninds), theta(:,k));
        [lambdaA, dlambdaA] = Combined_rateCalcs(allX(1,allinds), allX(2,allinds), theta(:,k));

        
        dlambdaA_lambdaA = zeros(6,6);
        for m = 1:6
            for n = m:6
                dlambdaA_lambdaA(m,n) = sum(dlambdaA(m,:).*dlambdaA(n,:).*lambdaA)*samplingInterval;
                if (m < n)
                    dlambdaA_lambdaA(n,m) = dlambdaA_lambdaA(m,n);
                end
                
            end
        end
        dll = sum(dlambdaT,2) - samplingInterval*sum(dlambdaA.*[lambdaA;lambdaA;lambdaA;lambdaA;lambdaA;lambdaA],2);
        w(:,:,j) = inv( inv(w(:,:,k) + q_rate*dt(k)) + dlambdaA_lambdaA);
        theta(:,j) = theta(:,k) + w(:,:,j)*dll;
        diagnostics.innovation(j) = nnz(turninds) - sum(lambdaA)*samplingInterval;
        diagnostics.innovation_scaled(j) = diagnostics.innovation(j)/sqrt(sum(lambdaA)*samplingInterval);
        na = nnz(allinds);
        diagnostics.innovation_scaled_alt(j) = (sum(1./(lambdaT*samplingInterval)) - na)/sqrt(na);
        diagnostics.dll(:,j) = dll;
        diagnostics.dT(:,j) = w(:,:,j)*dll;
    end
    
else % for one (1) white noise stimulus
    
    
    turnTime = turnTime(:)';
    turnX = turnX(:)';
    runTime = runTime(:)';
    runX = runX(:)';
    
    allTime = [turnTime runTime];
    allX = [turnX runX];
    
    if (nargin < 6)
        samplingInterval = median(diff(allTime));
    end
    if (nargin < 7)
        tx =min(allTime):samplingInterval:max(allTime);
    end
    dt = diff(binEdgesFromCenters(tx));
    
    [~,~,turnBin] = histcounts(turnTime,  binEdgesFromCenters(tx));
    [~,~,allBin] = histcounts(allTime,  binEdgesFromCenters(tx));
    isrun = true(size(allBin));
    isrun(1:length(turnTime)) = true;
    
    w = zeros(3,3,length(dt));
    if (~exist('initW', 'var'))
        w(:,:,1) = diag([1 1 1]);
    else
        w(:,:,1) = initW;
    end
    theta = zeros(3,length(dt));
    if (~exist('initTheta', 'var'))
        v0 = 1;
        u0 = 0;
        a0 = log(length(turnTime)) - log(length(allTime)*samplingInterval);
       % theta(:,1) = [v0;u0;a0]; % for rateCalcs
        theta(:,1) = [0; 0; log(length(turnTime)) - log(length(allTime)*samplingInterval)]; %for rate2Calcs
    else
        theta(:,1) = initTheta;
    end
    %theta = repmat([a0;v0;u0], [1 length(dt)]);
    
    if (size(q_rate,1) ~= 3 || size(q_rate,2) ~= 3)
        q_rate = q_rate(1)*eye(3);
    end
    
    for j = 1:length(dt)
        k = max(j-1,1);
        turninds = turnBin == j;
        allinds = allBin == j;
        [lambdaT, dlambdaT, ddlambdaT] = rate2Calcs(turnX(turninds), theta(:,k));
        [lambdaA, dlambdaA, ddlambdaA] = rate2Calcs(allX(allinds), theta(:,k));
        ddlambdaT = sum(ddlambdaT,2);
        ddlambdaA_lambdaA = sum(ddlambdaA.*([lambdaA;lambdaA;lambdaA]), 2)*samplingInterval;
        ddl = ddlambdaT + ddlambdaA_lambdaA;
        dlambdaA_lambdaA = zeros(3,3);
        for m = 1:3
            for n = m:3
                dlambdaA_lambdaA(m,n) = sum(dlambdaA(m,:).*dlambdaA(n,:).*lambdaA)*samplingInterval;
                if (m < n)
                    dlambdaA_lambdaA(n,m) = dlambdaA_lambdaA(m,n);
                end
                
            end
        end
        dll = sum(dlambdaT,2) - samplingInterval*sum(dlambdaA.*[lambdaA;lambdaA;lambdaA],2);
        w(:,:,j) = inv( inv(w(:,:,k) + q_rate*dt(k)) - [0 ddl(1) 0; ddl(1) ddl(2) 0; 0 0 0] + dlambdaA_lambdaA);
        theta(:,j) = theta(:,k) + w(:,:,j)*dll;
        
        diagnostics.innovation(j) = nnz(turninds) - sum(lambdaA)*samplingInterval;
        diagnostics.innovation_scaled(j) = diagnostics.innovation(j)/sqrt(sum(lambdaA)*samplingInterval);
        na = nnz(allinds);
        diagnostics.innovation_scaled_alt(j) = (sum(1./(lambdaT*samplingInterval)) - na)/sqrt(na);
        diagnostics.dll(:,j) = dll;
        diagnostics.dT(:,j) = w(:,:,j)*dll;
    end
   

    
end
end

function [lambda, dlambda, ddlambda] = Combined_rateCalcs(x1, x2, theta)
lambda = exp( theta(1)*x1.^2 + theta(2)*x2.^2 + theta(3)*x1.*x2 + theta(4)*x1 + theta(5)*x2 + theta(6) );
if (nargout > 1)
    dlambda = [x1.^2; x2.^2; x1.*x2; x1; x2; ones(size(x1))];
end
if (nargout > 2)
    ddlambda = zeros(size(x1));
end
end

function [lambda, dlambda, ddlambda] = rateCalcs(x, theta)
lambda = exp(theta(3) + 0.5 * x.^2 - 0.5*theta(1).*(x-theta(2)).^2);
if (nargout > 1)
    dlambda = [-0.5.*(x-theta(2)).^2; theta(1).*(x-theta(2)); ones(size(x))];
end
if (nargout > 2)
    ddlambda = [(x-theta(2)); -theta(1)*ones(size(x)); zeros(size(x))];
end
end

function [lambda, dlambda, ddlambda] = rate2Calcs(x, theta)
lambda = exp(theta(3) +  theta(1).* x.^2 + theta(2).*x );
if (nargout > 1)
    dlambda = [x.^2; x; ones(size(x))];
end
if (nargout > 2)
    ddlambda = [zeros(size(x)); zeros(size(x)); zeros(size(x))];
end
end