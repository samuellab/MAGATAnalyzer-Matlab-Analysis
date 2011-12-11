function addGlobalQuantity(expt,varargin)
%adds globalQuantity to list of quantities experiment is keeping track of
%function addGlobalQuantity(expt,varargin)
%
%expt.addGlobalQuantity(globalQuantity)
%expt.addGlobalQuantity(xfieldname, newfieldname, xdata, ydata)
%expt.addGlobalQuantity(xfieldname, newfieldname, xdata, ydata, method)
%
%adds globalQuantity to list of quantities experiment is keeping track of
%(or creates a new global quantity from input params), then assigns values
%to tracks
%
%outputs: none
%inputs:
%EXPT: a member of the Experiment class 
%GLOBALQUANTITY: a GlobalQuantity created and defined elsewhere (see
%GlobalQuantity)
%or
%XFIELDNAME: the name of the field on which the global quantity depends
%   for instance, if the temperature is varying with time, xfieldname would
%   be 'eti'
%   if light intensity is a function of space, xfieldname might be 'sloc'
%NEWFIELDNAME: the name of the new field created, e.g. 'lightIntensity'
%XDATA: the values of xfieldname over which the new field function is
%       defined.
%       for instance if the temperature was measured every second, xdata
%       would be 0:1:lastTime
%YDATA: the values of the new field at each point specified by xdata
%METHOD: a function of the form
%        yfield = METHOD (xfield, xData, yData); see GlobalQuanity for
%        details
if isempty(varargin)
    warning('GERSHOW:AGC01', 'You have to pass something to addGlobalQuantity');
    return;
end

if isa(varargin{1}, 'GlobalQuantity')
    gq = varargin{1};
else
    gq = GlobalQuantity;
    fieldlist = {'xField','fieldname','xData','yData','derivationMethod'};
    n = min(length(fieldlist), length(varargin));
    if (n < 4)
        warning('GERSHOW:AGC01', 'I need at least 4 arguments to create a new GQ');
        return;
    end
    for j = 1:n
        gq.(fieldlist{j}) = varargin{j};
    end
end

if (isempty(expt.globalQuantity))
    expt.globalQuantity = gq;
else
    ind = find(strcmpi(gq.fieldname, {expt.globalQuantity.fieldname}));
    if (isempty(ind))    
        expt.globalQuantity = [expt.globalQuantity gq];
    else
        expt.globalQuantity(ind) = gq;
    end
end

for j = 1:length(expt.track)
    gq.addQuantityToTrack(expt.track(j));
end
