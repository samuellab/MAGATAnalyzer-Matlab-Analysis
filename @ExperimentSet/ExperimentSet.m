classdef ExperimentSet < handle
    %a class for manipulating groups of experiments
    %
    
    properties
        expt; % array of experiments
        defaultTitle = 'untitled'; % the default title prepended to auto generated graphs
        autocorr_tau = 0;
    end
    
    methods
        qvec = gatherField(eset, fieldname, varargin)
        qvec = gatherSubField (eset, field, subfield, varargin)
        qvec = gatherFromSubField(eset, subfield, fieldname, varargin)
        varargout = executeTrackFunction(eset, func, varargin);
        result = evaluateTrackExpression(eset, expression);
        varargout = executeExperimentFunction(eset, func, varargin);
        [h, eb] = makeHistogram(eset, fieldname, fieldaxis, varargin);
        [h, eb] = makeSubFieldHistogram(eset, field, subfield, fieldaxis, varargin);
        [h, eb] = makeReorientationHistogram(eset, fieldname, fieldaxis, varargin);
        [h, eb] = makeHeadSwingAcceptanceHistogram(eset, fieldname, fieldaxis, varargin);
        [x,meany,standarderror,standarddeviation] = meanField2vsField1 (eset, field1name, field2name, field1axis, varargin);
        [x,meany,standarderror,standarddeviation] = meanField2vsField1_slidingwindow (eset, field1name, field2name, field1center, field1binsize, windowType, varargin);
        [x,meany,standarderror,standarddeviation] = meanSubField2vsSubField1 (eset, field1name, subfield1name, field2name, subfield2name, field1axis, varargin);
        [xc, np, tx, nt] = crosscorrelate (eset, fieldname1, fieldname2, varargin);
        [ac, np, tx, nt] = autocorrelate (eset, fieldname, varargin);
        tau = getAutocorrTau(eset, varargin);
        setAutocorrTau(eset, varargin);
        [h,eb] = makeReorientationHistogram_slidingWindow(eset, fieldname, fieldcenters, fieldwidth, windowType, varargin)


        
        varargout = indToTrack(eset, ind);
        toMatFiles(eset, fstub);
        addTimingByFrameDifference(eset, deltaT);
    end
    methods(Static)
        eset = fromFiles(varargin);
        eset = fromMatFiles(fstub, fileinds,segment);
        eset = loadTrimStitchAndSave(basedir, esetname, ecl, camcalinfo, varargin);
    end
    
end

