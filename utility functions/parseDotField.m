function value = parseDotField(obj, dotfieldstring)
%function value = parseDotField(obj, dotfieldstring)

fields = textscan(dotfieldstring,'%s','Delimiter','.');
value = getfield(obj,fields{1}{:});