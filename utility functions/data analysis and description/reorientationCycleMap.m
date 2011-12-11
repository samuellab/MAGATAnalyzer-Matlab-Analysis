function rowdata = reorientationCycleMap (object, cyclic_field_name)
%function rowdata = reorientationCycleMap (object, cyclic_field_name)
%
%useful for generating raster style plots
rowdata = [];

if (length(object) > 1)
    for j = 1:length(object)
        rowdata = [rowdata reorientationCycleMap(object(j), cyclic_field_name)];%#ok<*AGROW>
    end
    return;
end
            
if (isa (object, 'ExperimentSet'))
    rowdata = reorientationCycleMap(object.expt, cyclic_field_name);
    return;
end

if (isa (object, 'Experiment'))
    rowdata = reorientationCycleMap(object.track, cyclic_field_name);
    return;
end

if (~isa(object, 'Track'))
    warning ('rcm:type','reorientationCycleMap expects an experiment set, experiment, or track');
    return;
end

t = object;

et = t.getDerivedQuantity('eti');
cet = t.getDerivedQuantity(cyclic_field_name);

run_start = [t.run.startInd];
run_end = [t.run.endInd];
reo_start = [t.reorientation.startInd];
reo_end = [t.reorientation.endInd];

if (isa (t, 'MaggotTrack'))
    hs_start = [t.headSwing.startInd];
    hs_end = [t.headSwing.endInd];
    hsa = [t.headSwing.accepted];
    nhs = [t.reorientation.numHS];
end

startInd = find(cet >= 0, 1, 'first');
endInd = find(diff(cet(startInd:end)) < 0, 1, 'first');
firstrunind = find([t.isrun], 1, 'first');
lastrunind = find([t.isrun], 1, 'last');
j = 1;
while (~isempty(endInd))
    endInd = startInd + endInd - 1;
    
    rowdata(j).trackNum = t.trackNum; 
    rowdata(j).allvalid = cet(startInd) <= 2*t.dr.interpTime & endInd < length(et);
    rowdata(j).startTime = et(startInd) - cet(startInd);
    rowdata(j).endTime = et(endInd);
    rowdata(j).validStart = cet(startInd);
    rowdata(j).validEnd = cet(endInd);
    rowdata(j).startsRunning = logical(t.isrun(startInd));
    rowdata(j).run_start = cet(run_start(run_start >= startInd & run_start <= endInd));
    
    rowdata(j).run_end = cet(run_end(run_end >= startInd & run_end <= endInd));
    
    reoinds = find(reo_start >= startInd & reo_end <= endInd);
    firstreoind = find(reo_start < startInd & reo_end > startInd);
    lastreoind = find(reo_end > endInd & reo_start < endInd);
    
    rowdata(j).reo_start = cet(reo_start(reoinds));
    rowdata(j).reo_end = cet(reo_end(reoinds));
    
    if (~isempty(firstreoind))
        rowdata(j).reo_start = [cet(startInd) rowdata(j).reo_start];
        rowdata(j).reo_end = [cet(reo_end(firstreoind)) rowdata(j).reo_end];
    end
     if (~isempty(lastreoind))
        rowdata(j).reo_start = [rowdata(j).reo_start cet(reo_start(lastreoind))];
        rowdata(j).reo_end = [rowdata(j).reo_end cet(endInd)];
    end
    
%     rowdata(j).reo_start = cet(reo_start(reo_start >= startInd & reo_start <= endInd));
%     rowdata(j).reo_end = cet(reo_end(reo_end >= startInd & reo_end <= endInd));
%    
    if (startInd < firstrunind)
        rowdata(j).prerun_start = cet(startInd);
        rowdata(j).prerun_end = cet(min(firstrunind, endInd));
    else
        rowdata(j).prerun_start = [];
        rowdata(j).prerun_end = [];
    end
    if (lastrunind < endInd)
        rowdata(j).postrun_start = cet(max(lastrunind,startInd));
        rowdata(j).postrun_end = cet(endInd);
     else
        rowdata(j).postrun_start = [];
        rowdata(j).postrun_end = [];
    end

    if (isa (t, 'MaggotTrack'))
        rowdata(j).ahs_start = cet(hs_start(hs_start >= startInd & hs_start <= endInd & hsa));
        rowdata(j).rhs_start = cet(hs_start(hs_start >= startInd & hs_start <= endInd & ~hsa));
        rowdata(j).ahs_end = cet(hs_end(hs_end >= startInd & hs_end <= endInd & hsa));
        rowdata(j).rhs_end = cet(hs_end(hs_end >= startInd & hs_end <= endInd & ~hsa));
%         rowdata(j).reo_start_nhs = nhs(reo_start >= startInd & reo_start <= endInd);
%         rowdata(j).reo_end_nhs = nhs(reo_end >= startInd & reo_end <= endInd);
        rowdata(j).reo_start_nhs = nhs([firstreoind reoinds lastreoind]);
        rowdata(j).reo_end_nhs = rowdata(j).reo_start_nhs;
    end
%      if (~isempty(rowdata(j).reo_end) && ~isempty(rowdata(j).reo_start) &&...
%             rowdata(j).reo_end(1) < rowdata(j).reo_start(1) && rowdata(j).run_start(1) > rowdata(j).reo_end(1))
%         rowdata(j).reo_start = [cet(startInd) rowdata(j).reo_start];
%         rowdata(j).reo_start_nhs(1) = rowdata(j).reo_end_nhs(1);
%     end
%     if (~isempty(rowdata(j).reo_end) && ~isempty(rowdata(j).reo_start) &&...
%             rowdata(j).reo_start(end) > rowdata(j).reo_end(end) && rowdata(j).run_end(end) < rowdata(j).reo_start(end))
%         rowdata(j).reo_end = [rowdata(j).reo_end cet(endInd)];
%         rowdata(j).reo_end_nhs(end) = rowdata(j).reo_start_nhs(end);
%     end
    rowdata(j).endsRunning = logical(t.isrun(endInd));
    startInd = endInd + 1;
    endInd = find(diff(cet(startInd:end)) < 0, 1, 'first');
    j = j+1;
    if (isempty(endInd) && startInd < length(et))
        endInd = length(et) + 1 - startInd;
    end
end
% if (startInd < length(et))
%     endInd = length(et);
%     rowdata(j).trackNum = t.trackNum; 
%     rowdata(j).allvalid = false;
%     rowdata(j).startTime = et(startInd) - cet(startInd);
%     rowdata(j).endTime = et(endInd);
%     rowdata(j).validStart = cet(startInd);
%     rowdata(j).validEnd = cet(endInd);
%     rowdata(j).startsRunning = logical(t.isrun(startInd));
%     rowdata(j).run_start = cet(run_start(run_start >= startInd & run_start <= endInd));
%     rowdata(j).reo_start = cet(reo_start(reo_start >= startInd & reo_start <= endInd));
%     rowdata(j).run_end = cet(run_end(run_end >= startInd & run_end <= endInd));
%     rowdata(j).reo_end = cet(reo_end(reo_end >= startInd & reo_end <= endInd));
%     
%     if (isa (t, 'MaggotTrack'))
%         rowdata(j).ahs_start = cet(hs_start(hs_start >= startInd & hs_start <= endInd & hsa));
%         rowdata(j).rhs_start = cet(hs_start(hs_start >= startInd & hs_start <= endInd & ~hsa));
%         rowdata(j).ahs_end = cet(hs_end(hs_end >= startInd & hs_end <= endInd & hsa));
%         rowdata(j).rhs_end = cet(hs_end(hs_end >= startInd & hs_end <= endInd & ~hsa));
%         rowdata(j).reo_start_nhs = nhs(reo_start >= startInd & reo_start <= endInd);
%         rowdata(j).reo_end_nhs = nhs(reo_end >= startInd & reo_end <= endInd);
%     end
%     if (startInd < firstrunind)
%         rowdata(j).prerun_start = cet(startInd);
%         rowdata(j).prerun_end = cet(min(firstrunind, endInd));
%     else
%         rowdata(j).prerun_start = [];
%         rowdata(j).prerun_end = [];
%     end
%     if (lastrunind < endInd)
%         rowdata(j).postrun_start = cet(max(lastrunind,startInd));
%         rowdata(j).postrun_end = cet(endInd);
%      else
%         rowdata(j).postrun_start = [];
%         rowdata(j).postrun_end = [];
%     end
%     rowdata(j).endsRunning = logical(t.isrun(endInd));
% end