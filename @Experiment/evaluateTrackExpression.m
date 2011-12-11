function result = evaluateTrackExpression(expt, expression)
%evaluates the expression for each track in experiment
%function result = evaluateTrackExpression(expt, expression)
%
%evaluates the expression for each track in experiment
%expression should use 'track' as the name of the track
%and 'expt' as the name of the experiment
%
%examples:
%to get the number of runs in each track
%numruns = expt.evaluateTrackExpression('length(track.run)');
%
%to set all the segmentation options to a default value
%expt.evaluateTrackExpression('track.so = expt.so');
%
%result is returned as a matrix if possible, otherwise as a cell
%
%outputs:
%RESULT: optional; the result of the expression for each track
%       asking for an output from an expression that doesn't produce one
%       will result in errors
%inputs:
%EXPT: a member of the experiment class
%EXPRESSION: the expression to evaluate

for j = 1:length(expt.track)
    track = expt.track(j); %#ok<NASGU>
    if (nargout > 0)
        mycell{j} = eval(expression); %#ok<AGROW>
    else
        eval(expression);
    end
end

if (nargout > 0)
    try 
        result = cell2mat(mycell);
    catch me
        result = mycell;
    end
end