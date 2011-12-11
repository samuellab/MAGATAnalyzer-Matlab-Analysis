function [x,meany,standarderror,standarddeviation] = meanSubField2vsSubField1 (eset, field1name, subfield1name, field2name, subfield2name, field1axis, varargin)
%function [x,meany,standarderror,standarddeviation] = meanSubField2vsSubField1 (eset, field1name, subfield1name, field2name, subfield2name, field1axis, varargin)
%
%gathers two sub fields (with varargin arguments for both) and calls meanyvsx on them
%if no return arguments are requested, plots the result
%note that field1axis is bin edges [edge1,edge2) and x is the mean value
%of subfield1 in that bin
%
%passing 'inds', inds will take value only from inds in this
%fashion (same for ydata/field2name):
%xdata = eset.gatherSubField(field1name, subfield1name, varargin{:});
%xdata = xdata(inds)

inds = [];
varargin = assignApplicable(varargin);


xdata = eset.gatherSubField(field1name, subfield1name, varargin{:});
ydata = eset.gatherSubField(field2name, subfield2name, varargin{:});

if (~isempty(inds))
    xdata = xdata(inds); %xdata must be 1-D
    ydata = ydata(:,inds);
end
[x,meany,standarderror,standarddeviation] = meanyvsx(xdata, ydata, field1axis);

if (nargout == 0)
    errorbar (x,meany,standarderror); title (eset.defaultTitle);
    xlabel(subfield1name);
    ylabel(['<' subfield2name '>']);
    embiggen();
end