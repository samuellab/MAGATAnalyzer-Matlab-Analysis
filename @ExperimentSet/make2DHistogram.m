function h2 = make2DHistogram(eset, fieldname1, fieldaxis1, fieldname2, fieldaxis2, varargin)
%function h2 = make2DHistogram(eset, fieldname1, fieldaxis1, fieldname2, fieldaxis2, varargin)
%
%generates a histogram of all values in expt.track.getDerivedQuantity(fieldname,varargin{:})
%if no arguments are specified, generates a plot of that histogram
%optional args 'subfield', subfieldname
%   calls eset.gatherFromSubfield(subfieldname, fieldname1,varargin{:});
subfield = [];
varargin = assignApplicable(varargin);
if (isempty(subfield))
    h = makeIm (eset.gatherField(fieldname1, varargin{:}), eset.gatherField(fieldname2, varargin{:}), fieldaxis1, fieldaxis2);
else
    h = makeIm (eset.gatherFromSubField(subfield, fieldname1, varargin{:}), eset.gatherFromSubField(subfield,fieldname2, varargin{:}), fieldaxis1, fieldaxis2);
end
if (nargout == 0)
    pcolor (fieldaxis1, fieldaxis2, h); shading flat; colorbar vert
    
    xlabel(fieldname1); ylabel(fieldname2); embiggen();
else
    h2 = h;
end
