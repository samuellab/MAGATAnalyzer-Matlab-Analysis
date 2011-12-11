function timinfo = createTimFile (directory, file_stub, outputfile, varargin)
%function timinfo = createTimFile (directory, file_stub, outputfile, varargin)
%
%creates a fake timing file based on file save times
%
%directory - the directory containing all the image files, including the
%last '\'
%file_stub - the part of the file preceding the numeric portion
%outputfile - the file to which timing info should be saved (to create a
%fake .tim file) -- optional
%timinfo - if you ask for a return, it will give you the three columns
%(including a column of 0s) that will be saved in the .tim file
%
%optional parameter/value pairs
%'extension', '.ext' -- if the extension is not jpg, pass it (including the
%dot) 
%'verbose', true/false -- if you want extra details
%
%example:
%dirfn = '\\labnas1\share\David\Image Data\david\20090306\18-23_C25_1\';
%fstub = 'N2g15_';
%createTimFile
verbose = false;
extension = '.jpg';
varargin = assignApplicable(varargin);

if (verbose)
    ts = tic;
    disp ('getting directory information');
end

d = dir([directory file_stub '*' extension]);

if (verbose)
    toc(ts);
    disp ([num2str(length(d)) ' files found']);
end

[inds,I] = sort(sscanf([d.name], [file_stub '%d' extension]));
inds = inds';
pictim = [d(I).datenum]; %datenum is in days
pictim = reshape(pictim, size(inds));
pictim = (pictim-min(pictim))* 24 * 3600;
%stepinds = find(diff(pictim) > 0);
%sigma = 3 * max(diff(stepinds));
sigma = 10;
if (verbose)
    disp (['low passing with sigma = ' num2str(sigma)]);
end
mselapsed = 1000 * lowpass1D(pictim, sigma,'padType','mirror');

triggernum = zeros(size(inds));
tinfo = uint32(round([mselapsed;triggernum;pictim]));

existsAndDefault('outputfile', []);
if (~isempty(outputfile))
    if (verbose)
        disp (['saving tim info to ', outputfile]);
    end
    fid = fopen(outputfile, 'w');
    fprintf(fid, '%d\t%d\t%d\n', tinfo);
    fclose(fid);
end
if (nargout > 0)
    timinfo = tinfo;
end