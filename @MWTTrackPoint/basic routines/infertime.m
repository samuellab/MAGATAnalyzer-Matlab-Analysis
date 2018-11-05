function timeinsecs = infertime (fstub, ext)
%function timeinsecs = infertime (fstub, ext)
%
%infers time picture was taken from timestamp on files

d = dir ([fstub '*.' ext]);
t = [d.datenum];
t = sort(t);
t = t - t(1);
timeinsecs = t*3600*24;
%timeinsecs = lowpass1D(t*3600*24,10);

