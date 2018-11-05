function sdd = findSuppDataDir(exp_fname)
%function sdd = findSuppDataDir(exp_fname)
%
%finds supplemental data directory for gershow lab style processing

 [p,f] = fileparts(exp_fname);
 sdd = fixFileNameWin(fullfile(p, [f ' sup data dir']));
 if (isdir(sdd))
     return;
 end
 
 sdd = fixFileNameWin(fullfile(p, [f ' sup data']));
 if (isdir(sdd))
     return;
 end
 
 sdd = fixFileNameWin(fullfile(p, [shortenFileStub(f) ' sup data dir']));
 if (isdir(sdd))
     return;
 end
 
 sdd = fixFileNameWin(fullfile(p, [shortenFileStub(f) ' sup data']));
 if (isdir(sdd))
     return;
 end
 
 sdd = '';
 
