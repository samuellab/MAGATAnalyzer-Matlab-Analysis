function  str = getReport(run, varargin)
%generates a report about a run
%function  getReport(run, varargin)
%
%output: STR (a string, or cell of strings if multiple runs)
%input: run < Run
%optional arguments:
%   meanFields: default {'speed'}: fields over which to display a mean
%   value

str = ['Run: ', getReport@TrackPart(run, varargin{:})];
