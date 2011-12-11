function info = createOutputFileNames (info, inputbasedir, outputbasedir, imdirname, extdirname)
%function info = createOutputFileNames (info, inputbasedir, outputbasedir, imdirname, extdirname)
%
% INFO: structure created by parseDirectoryAndFileName
% INPUTBASEDIR (optional): the base directory where images are stored, this part of
% the path will be replaced with outputbasedir, assumed to be at the beginning
% of the file path
% OUTPUTBASEDIR (required if inptubasedir is defined): the base directory where extracted images are to be placed
% IMDIRNAME (optional): the name of the folder holding image data; this will be
% changed to extdirname
% EXTDIRNAME (required if imdirname is defined): the name of the folder
% holding extracted data
% pass empty brackets to bypass arguments
%
% example for style 'marc'
% image_stack is 
% \\labnas2\LarvalCO2\Image Data\ethyl acetate 2000X diluted Spatial\50mL air through 2000X diluted ethyl acetate\CS\20110215\1657\CS_stack.mmf
% info = parseDirectoryAndFileName(image_stack) results in:
% info.descriptiveStub = 'CS'
% info.timeStamp = '1657'
% info.date = '20110215'
% info.genotype = 'CS'
% info.experimentType = '50mL air through 2000X diluted ethyl acetate'
% info.basepath = '\\labnas2\LarvalCO2\Image Data\ethyl acetate 2000X diluted Spatial\'
%
% inputbasedir = '\\labnas2\LarvalCO2\Image Data'
% outputbasedir = 'E:\LarvalCO2\Extracted'
% imdirname = []
% extdirname = []
% or
% inputbasedir = '\\labnas2\'
% outputbasedir = 'E:\'
% imdirname = 'Image Data'
% extdirname = 'Extracted'
%
% both result in extracted bin file:
%'E:\LarvalCO2\Extracted\ethyl acetate 2000X diluted Spatial\50mL air through 2000X diluted ethyl acetate\CS\20110215_1657_CS_tracks.bin'
anythingBueller = false;
anythingBueller = anythingBueller | (existsAndDefault ('inputbasedir', []) & existsAndDefault ('outputbasedir', []));
anythingBueller = anythingBueller | (existsAndDefault ('imdirname', []) & existsAndDefault ('extdirname', []));
if (~anythingBueller)
    disp ('type help createOutputFileNames for more info');
    error ('at least one pair of inputbasedir/outputbasedir or imdirname/extdirname must be defined');
end

switch (lower(info.style))
    case 'marc'
        info.outputbasepath = swapOutDirectories (info.basepath, inputbasedir, outputbasedir, imdirname, extdirname);
        info.outputstub = fullfile(info.outputbasepath, info.experimentType, info.genotype, [info.date '_' info.timeStamp '_' info.descriptiveStub '_tracks']);
        info.outputBinFile = [info.outputstub '.bin'];
        info.headerFile = [info.outputstub '_header.txt'];
        info.logFile = [info.outputstub '_extraction_log.txt'];
    otherwise
        error (['style: ' info.style ' not defined.  Only ''marc'' is currently implemented']);
end
        
function outputpath = swapOutDirectories (basepath, inputbasedir, outputbasedir, imdirname, extdirname)
    if (~isempty(inputbasedir))
        [common, residual] = getCommonPath(basepath, inputbasedir);
        if (isempty(common))
            warning (['basepath: ' basepath ' and inputbasedir ' inputbasedir ' have nothing in common']);
        else
            if (~isSameFile(common, inputbasedir))
                basepath
                inputbasedir
                outputbasedir
                imdirname
                extdirname
                [c, r] = getCommonPath(inputbasedir, common);
                warning ([c ' is common to base path and inputbasedir, but ' r ' is found only in inputbasedir']);
            end
        end
        
        outputpath = fullfile(outputbasedir, residual);
    else
        outputpath = basepath;
    end
    
    if (~isempty(imdirname))
        c = pathToCellOfStrings(outputpath);
        ind = find(strcmpi(imdirname, c));
        if (isempty(ind))
            warning ([imdirname ' is not parth of base path: ' basepath]);
        else
            if (length(ind) > 1)
                warning ([imdirname ' appears multiple times in base path: ' basepath ' (replacing first instance only)']);
                ind = ind(1);
            end
            c{ind} = extdirname;
        end
        outputpath = fullfile(c{:});
    end
    
    if (isSameFile(basepath, outputpath))
        warning ('basepath and outputpath are identical');
    end
