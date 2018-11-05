function fs = shortenFileStub (f)
%function fs = shortenFileStub (f)
%
%shortens a file stub from an experiment file name by removing all 
%text between the first and last underscores ( _ )
%if there are not two underscores or if the resulting stub is too long
%takes at most the first 12 and last 12 characters from the file name
%

[~,f,~] = fileparts(f);

t = regexp(f, '_', 'split');
t1 = t{1};
t2 = t{end};
t1 = t1(1:min(length(t1),12));
t2 = t2(max(1, length(t2)-11):end);
fs = [t1 '_' t2];
