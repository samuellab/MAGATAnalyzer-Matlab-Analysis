function [sigma trialValues] = stdOfFunVal(modelParams, modelLogLikelihoodHessian, funOfModelParams, varargin)
%function [minval maxval] = confidenceRange(modelParams, modelLogLikelihoodHessian, funOfModelParams, varargin)
%parameter value pairs:
%numtrials, 1E6 : number of points to calculate for taking the minimum and
%maximum
%
%modelParams: the parameters that maximize the likelihood of the model: (k x 1)
%modelLogLikelihoodHessian: the second derivative matrix at the maximum (k x k)
%funOfModelParams: should take a kxN vector of inputs and return a 1xN
%output

numtrials = 1E6;

varargin = assignApplicable(varargin);
if (size(modelParams, 1) == 1)
    modelParams = modelParams';
end
k = length(modelParams);

x = randn(k, numtrials);
if (any (~isfinite(modelLogLikelihoodHessian)))
    disp ('model hessian has non finite values')
    sigma = NaN;
    trialValues = [];
    return;
end

[V,D] = eig(modelLogLikelihoodHessian);
if (any(D > 0))
    disp ('model hessian has non-negative eigenvalues')
    sigma = NaN;
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
sigma = std(trialValues);



end

