function result = evaluateTrackExpression(eset, expression)
% evaluates a given expression at every track
% function result = evaluateTrackExpression(eset, expression)
%
% concatinates the results of eset.expt(j).evaluateTrackExpression(eset,
% expression) (if any)
% attempts to convert the resulting cell array to a matrix 
%
% details from Experiment/evaluateTrackExpression
% evaluates the expression for each track in experiment
% expression should use 'track' as the name of the track
% and 'expt' as the name of the experiment
%
% examples:
% to get the number of runs in each track
% numruns = expt.evaluateTrackExpression('length(track.run)');
%
% to set all the segmentation options to a default value
% expt.evaluateTrackExpression('track.so = expt.so');
%
% result is returned as a matrix if possible, otherwise as a cell
%
% outputs:
% RESULT: optional; the result of the expression for each track
%       asking for an output from an expression that doesn't produce one
%       will result in errors
% inputs:
% ESET: a member of the ExperimentSet class
% EXPRESSION: the expression to evaluate

for j = 1:length(eset.expt)
    if (nargout > 0)
        mycell{j} = eset.expt(j).evaluateTrackExpression(expression);
    else
        eset.expt(j).evaluateTrackExpression(expression);
    end
end


if (nargout > 0)
    try 
        result = cell2mat(mycell);
    catch me
        result = mycell;
    end
end


