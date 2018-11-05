function [minval maxval trialValues] = confidenceRange(modelParams, modelLogLikelihoodHessian, funOfModelParams, confidenceLevel, varargin)
%function [minval maxval] = confidenceRange(modelParams, modelLogLikelihoodHessian, funOfModelParams, confidenceLevel, varargin)
%parameter value pairs:
%numtrials, 1E6 : number of points to calculate for taking the minimum and
%maximum
%
%modelParams: the parameters that maximize the likelihood of the model: (k x 1)
%modelLogLikelihoodHessian: the second derivative matrix at the maximum (k x k)
%funOfModelParams: should take a kxN vector of inputs and return a 1xN
%output
if (length(confidenceLevel) > 1)
    for j = 1:length(confidenceLevel)
        [minvt maxvt tvt] = confidenceRange(modelParams, modelLogLikelihoodHessian, funOfModelParams, confidenceLevel(j), varargin{:});
        minval(j) = minvt;
        maxval(j) = maxvt;
        if (nargout > 2)
            trialValues(j,:) = tvt;
        end
    end
    return;
end

numtrials = 1E6;

varargin = assignApplicable(varargin);
if (size(modelParams, 1) == 1)
    modelParams = modelParams';
end
k = length(modelParams);
%{
prefactor = 1;
while true
    prefactor = prefactor*1.1;
    x = randn(k, ceil(prefactor*numtrials/confidenceLevel));
    maxrsq = chi2inv(confidenceLevel, k);
    x = x(:,sum(x.^2,1) < maxrsq);
    if (size(x,2) >= numtrials)
        break;
    end
end
x = x(:,1:numtrials);
%}
x = randomFromKBall(k, numtrials);
if (any (~isfinite(modelLogLikelihoodHessian)))
    disp ('model hessian has non finite values')
    minval = NaN;
    maxval = NaN;
    trialValues = [];
    return;
end

[V,D] = eig(modelLogLikelihoodHessian);
if (any(D > 0))
    disp ('model hessian has non-negative eigenvalues')
    minval = NaN;
    maxval = NaN;
    trialValues = [];
    return;
end

trialPoints = repmat(modelParams,1,size(x,2));
s = 1./sqrt(-D(logical(eye(size(D,1)))));
for j = 1:length(modelParams)
    trialPoints = trialPoints+s(j)*V(:,j)*x(j,:);
end
%whos trialPoints
trialValues = funOfModelParams(trialPoints);
%size(trialValues)
dim = find(size(trialValues) == length(trialPoints));
if (isempty(dim))
    dim = 1;
end
minval = min(trialValues,[],dim);
maxval = max(trialValues,[],dim);



end

