function pt = reloadPoint (expt, pt)
%loads a new copy of point from the position of pt in the bin file
%function pt = reloadPoint (expt, pt)
%
%outputs: 
%PT, a TRACKPOINT of same subclass as PT
%inputs:
%EXPT: a member of the Experiment class
%PT: a member of the TrackPoint class or a vector of TrackPoints

if (expt.fid == 0 || ftell(expt.fid) < 0)
    expt.openDataFile();
end


for j = 1:length(pt)
    fseek(expt.fid, pt(j).locInFile, -1);
    pt2(j) = pt(j).fromFile(expt.fid, true, true, expt.camcalinfo); %#ok<AGROW>
end
pt = pt2;