function f = fixFileNameWin(f)
%function f = fixFileNameWin(f)
%
%if is windows and is longer than 260 characters, generates short path (8
%char) directory names using windows API

MAX_PATH = 260;
if (~ispc || length(f) < MAX_PATH)
    return;
end

f = getshortpath(f);
