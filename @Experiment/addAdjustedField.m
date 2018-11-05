function addAdjustedField(expt, field1, adjfield,varargin) 
%function addAdjustedField(expt, field1, adjfield,varargin)
%
% EXPT < Experiment
% FIELD1 < 1-dimensional field (e.g. temperature)
% ADJFIELD < 1-dimensional field (e.g. speed)
%
% calculates mean adjfield vs. field1, then creates adjusted field to take out linear contribution
% say field1 = T; adjfield = s. 
% approximate <s(T)> = <s> + a(T-T0), where <s> is mean over all
% temperatures
% s_adj = s*F(T): F(T) = 1/(1+a/<s>*(T-T0));

x = expt.gatherField(field1);
y = expt.gatherField(adjfield);

inds = isfinite(x) & isfinite(y);
x = x(inds); y = y(inds);
x = [x;x]; x(2,:) = 1;
A = y*pinv(x);
my = mean(y);

xaxis = linspace(min(x(1,:)),max(x(1,:)), 100);
adj_vs_x = my./(A(1)*xaxis + A(2));

gq = GlobalQuantity;
gq.xField = {field1, adjfield};
gq.fieldname = [field1 '_adjusted_' adjfield];
gq.xData = xaxis;
gq.yData = adj_vs_x;
gq.derivationMethod = @GlobalQuantity.oneDinterpolationAndMultiplication;

expt.addGlobalQuantity(gq);