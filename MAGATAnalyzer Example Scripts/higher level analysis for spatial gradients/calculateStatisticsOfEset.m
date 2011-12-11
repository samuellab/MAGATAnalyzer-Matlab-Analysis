function eset_statistics = calculateStatisticsOfEset(eset, varargin )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if (length(eset) > 1)
    for j = 1:length(eset)
        es(j) = calculateStatisticsOfEset(eset(j), varargin{:}); %#ok<AGROW>
    end
    eset_statistics = es;
    return;
end

timerange = [];
validname = [];
validoperation = func2str(@(x) logical(setNonFiniteToZero(x)));
directions = [0 180 -90 90];
varargin = assignApplicable(varargin); %#ok<NASGU>

if (~isempty(timerange))
    validname = 'eti';
    validoperation =  ['@(x) x >= ' num2str(min(timerange)) ' & x <= ' num2str(min(timerange)) ];
end
es.validname = validname;
es.validoperation = validoperation;
es.directions = directions;

if (~isempty(validname) && ischar(validoperation))
    validoperation = str2func(validoperation);
end

es.directions = directions;


es.runStartTheta = eset.gatherSubField ('run', 'startTheta');
es.runEndTheta = eset.gatherSubField ('run', 'endTheta');
es.runMeanTheta = eset.gatherSubField('run' ,'meanTheta');

es.reoPrevDir = eset.gatherSubField ('reorientation', 'prevDir');
es.reoNextDir = eset.gatherSubField ('reorientation', 'nextDir');
es.reoNumHS = eset.gatherSubField ('reorientation', 'numHS');

es.headSwingAccepted = eset.gatherSubField ('headSwing', 'accepted');
es.headSwingPrevDir = eset.gatherSubField ('headSwing', 'prevDir');
es.firstHeadSwingPrevDir = eset.gatherSubField ('firsths', 'prevDir');
es.firstHeadSwingAccepted = eset.gatherSubField ('firsths', 'accepted');
es.lastHeadSwingAccepted = eset.gatherSubField ('lasths', 'accepted');
es.headSwingValid = eset.gatherSubField ('headSwing', 'valid');
es.firstHeadSwingValid = eset.gatherSubField ('firsths', 'valid');
es.lastHeadSwingValid = eset.gatherSubField ('lasths', 'valid');

if (~isempty(validname))
    es.runvfield_start = eset.gatherFromSubField ('run', validname, 'position', 'start');
    es.runvfield_end = eset.gatherFromSubField ('run', validname, 'position', 'end');
    es.reovfield = eset.gatherFromSubField ('reorientation', validname, 'position', 'start');
    es.hsvfield = eset.gatherFromSubField ('headSwing', validname, 'position', 'start');
    es.firsthsvfield = eset.gatherFromSubField ('firsths', validname, 'position', 'start');
    es.lasthsvfield = eset.gatherFromSubField ('lasths', validname, 'position', 'start');
    
    
    es.runStartValid = validoperation(es.runvfield_start);
    es.runEndValid = validoperation(es.runvfield_end);
    es.reoValid = validoperation(es.reovfield);
    es.headSwingValid = es.headSwingValid & validoperation(es.hsvfield);
    es.firstHeadSwingValid = es.firstHeadSwingValid & validoperation(es.firsthsvfield);
    es.lastHeadSwingValid = es.lastHeadSwingValid & validoperation(es.lasthsvfield);
else
    es.runvfield_start = [];
    es.runvfield_end = [];
    es.runStartValid = true(size(es.runStartTheta));
    es.runEndValid = true(size(es.runEndTheta));
    es.reoValid = true(size(es.reoPrevDir));
  
end

eti = eset.gatherField('eti');
if (isempty(timerange))
    timerange = [min(eti)-1 max(eti)+1];
end
es.timerange = timerange;

it = eset.gatherSubField('dr','interpTime');
dt = median(it);
if (any(it ~= dt))
    disp (['warning:  eset does not have homogenous interpolation times, instead range from ' num2str(min(it)) ' to ' num2str(max(it))]);
end

es.maxNumAnimals = zeros(size(eset.expt));
es.maxNumAnimalsValid = es.maxNumAnimals;
for j = 1:length(eset.expt)
    et = eset.expt(j).gatherField('eti');
    if (~isempty(validname))
        valid = validoperation(eset.expt(j).gatherField(validname));
    else
        valid = true(size(et));
    end

    es.maxNumAnimals(j) = max(histc(et, min(et):1:max(et)))*dt;
    es.maxNumAnimalsValid(j) = max(histc(et(valid), min(et):1:max(et)))*dt;
end
es.numExpts = length(eset.expt);
es.numAnimals = sum(es.maxNumAnimalsValid);
eti = eset.expt(j).gatherField('eti');
if (~isempty(validname))
    valid = validoperation(eset.gatherField(validname));
else
    valid = true(size(eti));
end
es.animalTime = nnz(valid)*dt;

for j = 1:length(directions)
    
    es.numRunsFromDirection(j) = nnz(cos(es.runStartTheta - deg2rad(directions(j))) > 1/sqrt(2) & es.runStartValid);
    es.numRunsInDirection(j) = nnz(cos(es.runMeanTheta - deg2rad(directions(j))) > 1/sqrt(2) & es.runStartValid);
    es.numReosFromDirection(j) = nnz(cos(es.reoPrevDir - deg2rad(directions(j))) > 1/sqrt(2) & es.reoValid);
    es.numReosWithHSFromDirection(j) = nnz(cos(es.reoPrevDir - deg2rad(directions(j))) > 1/sqrt(2) & es.reoValid & es.reoNumHS > 0);
    es.numHSFromDirection(j) = nnz(cos(es.headSwingPrevDir - deg2rad(directions(j))) > 1/sqrt(2) & es.headSwingValid);
end
if (~isempty(es.validoperation) && ~ischar(es.validoperation))
    es.validoperation = func2str(es.validoperation);
end
eset_statistics = es;

