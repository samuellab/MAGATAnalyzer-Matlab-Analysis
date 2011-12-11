function YamlStruct = ReadYaml(yaml_file)
% This function reads Yaml file into struct
% Example
% >> yaml_file = 'EnaspolMain.yaml';
% >> YamlStruct = ReadYaml(yaml_file)
%
%  %======================================================================
%{
		Copyright (c) 2011
		This program is a result of a joined cooperation of Energocentrum
		PLUS, s.r.o. and Czech Technical University (CTU) in Prague.
        The program is maintained by Energocentrum PLUS, s.r.o. and
        licensed under the terms of MIT license. Full text of the license
        is included in the program release.
		
        Author(s):
		Jiri Cigler, Dept. of Control Engineering, CTU Prague 
		Jan  Siroky, Energocentrum PLUS s.r.o.
		
        Implementation and Revisions:

        Auth  Date        Description of change
        ----  ---------   -------------------------------------------------
        jc    01-Mar-11   First implementation
        jc    02-Mar-11   .jar package initialization moved to external fun
        jc    18-Mar-11   Warning added when imported file not found
        jc    07-Jun-11   Ability to merge structures
        jc    09-Jun-11   Merging of structures recursively
%}
%======================================================================

InitYaml();

import('org.yaml.snakeyaml.Yaml');

yamlreader = Yaml();

Data = ReplaceImportByStruct(yaml_file);

work_folder=fileparts(yaml_file);
if isfield(Data,'import')
    yml = '';
    for i=1:numel(Data.import)
        fToImport=[work_folder filesep Data.import{i}];
        
        try            
            s = ReplaceImportByStruct(fToImport);
            Data = mergeStructs(Data,s);
        catch
            warning('YAMLMatlab:FileNotFoundException','YAMLMatlab: File %s not found',fToImport);
        end        
    end
    %jymlobj = yamlreader.load(yml);
    %YamlStruct = Hash2Struct(jymlobj);
    YamlStruct = rmfield(Data,'import');
else
    YamlStruct =Data;
end
    function s = ReplaceImportByStruct(fname)
        yml_file = fileread(fname);
        ymlobj = yamlreader.load(yml_file);

        s = Hash2Struct(ymlobj);
    end % end of ReplaceImportByStruct

end % end of function


function res = mergeStructs(x,y)

    
if isstruct(x) && isstruct(y)
    res = x;
    names = fieldnames(y);
    for fnum = 1:numel(names)
        if isfield(x,names{fnum})
            res.(names{fnum}) = mergeStructs(x.(names{fnum}),y.(names{fnum}));
        else
            res.(names{fnum}) = y.(names{fnum});
        end
    end
elseif isstruct(x) && iscell(y)
    found = 0;
    for i=1:numel(y)
        if isequal(x,y{i})
            found =1;
            break;
        end
    end
    if not(found)
        res = [x,y];
    else
        res = y;
    end
elseif  iscell(x) && isstruct(y)
    found = 0;
    for i=1:numel(x)
        if isequal(x{i},y)
            found =1;
            break;
        end
    end
    if not(found)
        res = [x,y];
    else
        res = y;
    end
elseif iscell(x) && iscell(y)
    
    for i=1:min( numel(x), numel(y))
        res{i} = mergeStructs(x{i},y{i});
    end
    
    if numel(x) > i
        res(i+1:numel(x)) = x(i+1:end);
    else
        res(i+1:numel(y)) = y(i+1:end);
    end
else
    res = y;
end
end
