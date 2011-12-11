function str = WriteYamlToString(matlab_struct)
%function str = WriteYamlToString(matlab_struct)
% This function writers struct into Yaml string
% Example
% >> yaml_file_old = 'EnaspolMain.yaml';
% >> yaml_file_new = 'EnaspolMain.yaml';
% >> YamlStruct = ReadYaml(yaml_file_old)
% >> WriteYaml(yaml_file_new,YamlStruct)

%======================================================================
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
        ----  ---------   --------------------------------------------------
        jc    01-Mar-11   First implementation
        jc    02-Mar-11   .jar package initialization moved to external fun

%}
%======================================================================

InitYaml();

import('org.yaml.snakeyaml.Yaml');

yaml = Yaml();

Data = Struct2Hashmap(matlab_struct);

str = char(yaml.dump(Data));

end % end of function









