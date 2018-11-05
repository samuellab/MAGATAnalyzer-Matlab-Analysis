function sa = commandLineToJavaStringArray(cl)
%function sa = commandLineToJavaStringArray(cl)
%
%converts a command line into a java string array

sa = strCell2JavaStrArray(regexp(cl, '\s+', 'split'));
