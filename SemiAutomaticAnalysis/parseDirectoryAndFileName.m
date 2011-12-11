function info = parseDirectoryAndFileName (fullfilename, style, varargin)
%function info = parseDirectoryAndFileName (fullfilename, style, varargin)
%
% used by semi-automatic analysis
%
% parses directory structure and file name to get information about the
% experimental setup and genotype of animal in the experiment
%
% this will be used to create output files during extraction
%
% style determines the behavior of this function;  currently only 'marc' is
% implemented
%
% example for style 'marc'
% image stack is 
% \\labnas2\LarvalCO2\Image Data\ethyl acetate 2000X diluted Spatial\50mL
% air through 2000X diluted ethyl acetate\CS\20110215\1657\CS_stack.mmf
% info.fstub = \\labnas2\LarvalCO2\Image Data\ethyl acetate 2000X diluted Spatial\50mL
% air through 2000X diluted ethyl acetate\CS\20110215\1657\CS_stack
% info.extension = '.mmf';
% info.descriptiveStub = 'CS'
% info.timeStamp = '1657'
% info.date = '20110215'
% info.genotype = 'CS'
% info.experimentType = '50mL air through 2000X diluted ethyl acetate'
% info.basepath = '\\labnas2\LarvalCO2\Image Data\ethyl acetate 2000X diluted Spatial\'

existsAndDefault('style', 'marc');
info.originalName = fullfilename;
info.style = style;
switch (lower(style))
    case ('marc')
        %directory structure is .../experiment type/genotype/YYYYMMDD/HHMM/genotypeOrSomethingDescriptive_stack.mmf

        [pathstr, name, info.extension] = fileparts(fullfilename);
        info.fstub = fullfile(pathstr, name);
        ind = strfind(name, '_stack');
        if (~isempty(ind))
            info.descriptiveStub = name(1:(ind-1));
        else
            info.descriptiveStub = name;
        end
        
        [pathstr, info.timeStamp] = fileparts(pathstr);
        [pathstr, info.date] = fileparts(pathstr);
%         [pathstr, info.month] = fileparts(pathstr);
%         [pathstr, info.year] = fileparts(pathstr);
        [pathstr, info.genotype] = fileparts(pathstr);
        [pathstr, info.experimentType] = fileparts(pathstr);
        info.basepath = pathstr;
    case ('janelia')
        %directory structure is
        %tracker/protocol/gal4(effector)/2011/09/27/1415/datafiles
        [pathstr, name, info.extension] = fileparts(fullfilename);
        info.fstub = fullfile(pathstr, name);
        info.descriptiveStube = name;
        [pathstr, info.timeStamp] = fileparts(pathstr);
        [pathstr, info.day] = fileparts(pathstr);
        [pathstr, info.month] = fileparts(pathstr);
        
        
        
    otherwise
        error (['style ' style ' not implemented yet']);
end
        
