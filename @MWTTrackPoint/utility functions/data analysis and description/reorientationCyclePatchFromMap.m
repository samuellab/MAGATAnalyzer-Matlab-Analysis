function [invalidpatches, runpatches, reopatches, ahspatches, rhspatches, prerunpatches, postrunpatches] = reorientationCyclePatchFromMap(rowdata, varargin)
%function [invalidpatches, runpatches, reopatches, ahspatches, rhspatches, prerunpatches, postrunpatches] = reorientationCyclePatchFromMap(rowdata, varargin)
%mintime = -Inf;
%maxtime = Inf;
%minhs = 1;
%maxhs = Inf;
%    include reorientation if nhs in [minhs,maxhs]
%includepartial = true;
%yoffset = 0;

%
% invalidcolor = [0.1 0.1 0.1];
% runcolor = [1 1 1];
% reocolor = [0.4 0.4 0.4];
% ahscolor = [0 1 0];
% rhscolor = [1 0 0];
mintime = -Inf;
maxtime = Inf;
minhs = 1;
maxhs = Inf;
includepartial = true;
yoffset = 0;

runpatches = [];
reopatches = [];
ahspatches = [];
rhspatches = [];
prerunpatches = [];
postrunpatches = [];

if (~isempty(varargin) && isstruct(varargin{1}))
    assignFromStruct(varargin{1});
    varargin = varargin{2:end};
else
    varargin = assignApplicable(varargin);
end

% colornames = {'invalidcolor', 'runcolor', 'reocolor', 'ahscolor', 'rhscolor'};
% for j = 1:length(colornames)
%     if (ischar(eval(colornames{j})))
%         eval([colornames{j} ' = char2rgb(' colornames{j} ');']);
%     end
% end

valid = false(size(rowdata));
reovalid = valid;
ahsvalid = valid;
rhsvalid = valid;
prerunvalid = valid;
postrunvalid = valid;
maxx = max([rowdata.validEnd]);
minx = min([rowdata.validStart]);
y = yoffset - 1;
for j = 1:length(rowdata)   
    if (~includepartial && ~rowdata(j).allvalid)
        continue;
    end
    if (rowdata(j).startTime < mintime || rowdata(j).endTime > maxtime)
        continue;
    end
    if (isempty(rowdata(j).run_start))
        if (rowdata(j).startsRunning)
            if (rowdata(j).allvalid)
                rs = minx;
            else
                rs = rowdata(j).validStart;
            end
        else
            continue;
        end
    else
        rs = rowdata(j).run_start;
    end
    if (isempty(rowdata(j).run_end))
        if (rowdata(j).endsRunning)
            if (rowdata(j).allvalid)
                re = maxx;
            else
                re = rowdata(j).validEnd;
            end
        else
            continue;
        end
    else
        re = rowdata(j).run_end;
    end

    
    valid(j) = true;
    y = y+1;
    
    if (~rowdata(j).allvalid)
        invalidpatches(j).Vertices = [minx y; rowdata(j).validStart y; rowdata(j).validStart y+1; minx y+1; rowdata(j).validEnd y; maxx y; maxx y+1; rowdata(j).validEnd y+1];
        invalidpatches(j).Faces = [1 2 3 4; 5 6 7 8];
    end
    
    
    if (re(1) < rs(1))
        if(~rowdata(j).startsRunning)
            j
            warning('run start/end alignment error');
        else
            rs = [minx rs];
        end
    end
    if (re(end) < rs(end))
        if(~rowdata(j).endsRunning)
            j
            warning('run start/end alignment error');
        else
            re = [re maxx];
        end
    end
    yverts = y*ones(size(rs))';
    if (length(rs) ~= length(re))
        j

        warning('run start length ~= run end length');
    else
        runpatches(j).Vertices = [rs' yverts; rs' yverts+1; re' yverts+1; re' yverts];
        N = length(rs);
        fv = (1:N)';
        runpatches(j).Faces = [fv fv+N fv+2*N fv+3*N];
    end
    rs = rowdata(j).reo_start(rowdata(j).reo_start_nhs >= minhs & rowdata(j).reo_start_nhs <= maxhs);
    re = rowdata(j).reo_end(rowdata(j).reo_end_nhs >= minhs & rowdata(j).reo_end_nhs <= maxhs);
    if (~isempty(rs) && ~isempty(re))
        reovalid(j) = true;
        if (re(1) < rs(1))
            if(rowdata(j).startsRunning)
                re = re(2:end);
                j
                warning('reo start/end alignment error');  %#ok<*WNTAG>
            end
            rs = [rowdata(j).validStart rs];
            
        end
        if (re(end) < rs(end))
            if(rowdata(j).endsRunning)
                rs(end)
                re(end)
                rs = rs(1:(end-1));
                j
                warning('reo start/end alignment error');
            end
            re = [re rowdata(j).validEnd];
        end
        yverts = y*ones(size(rs))';
        if (length(rs) ~= length(re))
            warning('reo start length ~= reo end length');
        end
        reopatches(j).Vertices = [rs' yverts; rs' yverts+1; re' yverts+1; re' yverts]; %#ok<*AGROW>
        N = length(rs);
        fv = (1:N)';
        reopatches(j).Faces = [fv fv+N fv+2*N fv+3*N];
    end
    rs = rowdata(j).ahs_start;
    re = rowdata(j).ahs_end;
    if (~isempty(rs) && ~isempty(re))
        
        if (re(1) < rs(1))
            if(rowdata(j).startsRunning)
%                 j
%                 warning('hs start/end alignment error');
            end
            rs = [rowdata(j).validStart rs];
            
        end
        if (re(end) < rs(end))
            if(rowdata(j).endsRunning)
%                 j
%                 warning('ahs start/end alignment error');
            end
            re = [re rowdata(j).validEnd];
        end
        yverts = y*ones(size(rs))';
        if (length(rs) ~= length(re))
%             j
%             warning('ahs start length ~= ahs end length');
        else
            ahspatches(j).Vertices = [rs' yverts; rs' yverts+1; re' yverts+1; re' yverts];
            N = length(rs);
            fv = (1:N)';
            ahspatches(j).Faces = [fv fv+N fv+2*N fv+3*N];
            ahsvalid(j) = true;
        end
    end
    rs = rowdata(j).rhs_start;
    re = rowdata(j).rhs_end;
    if (~isempty(rs) && ~isempty(re))
        
        if (re(1) < rs(1))
            if(rowdata(j).startsRunning)
%                 j
%                 warning('rhs start/end alignment error');
            end
            rs = [rowdata(j).validStart rs];
            
        end
        if (re(end) < rs(end))
            if(rowdata(j).endsRunning)
%                 j
%                 warning('rhs start/end alignment error');
            end
            re = [re rowdata(j).validEnd];
        end
        yverts = y*ones(size(rs))';
        if (length(rs) ~= length(re))
%             j
%             warning('rhs start length ~= rhs end length');
            
        else
            rhspatches(j).Vertices = [rs' yverts; rs' yverts+1; re' yverts+1; re' yverts];
            N = length(rs);
            fv = (1:N)';
            rhspatches(j).Faces = [fv fv+N fv+2*N fv+3*N];
            rhsvalid(j) = true;
        end
        rs = rowdata(j).prerun_start;
        re = rowdata(j).prerun_end;
        yverts = y*ones(size(rs))';
        if (~isempty(rs) && ~isempty(re))
            prerunpatches(j).Vertices = [rs' yverts; rs' yverts+1; re' yverts+1; re' yverts];
            N = length(rs);
            fv = (1:N)';
            prerunpatches(j).Faces = [fv fv+N fv+2*N fv+3*N];
            prerunvalid(j) = true;
        end
        rs = rowdata(j).postrun_start;
        re = rowdata(j).postrun_end;
        yverts = y*ones(size(rs))';
        if (~isempty(rs) && ~isempty(re))
            postrunpatches(j).Vertices = [rs' yverts; rs' yverts+1; re' yverts+1; re' yverts];
            N = length(rs);
            fv = (1:N)';
            postrunpatches(j).Faces = [fv fv+N fv+2*N fv+3*N];
            postrunvalid(j) = true;
        end
    end
end
if (~exist('invalidpatches', 'var'))
    invalidpatches = [];
else
    invalidpatches = invalidpatches(valid & ~[rowdata.allvalid]);
end
if (~isempty(runpatches)), runpatches = runpatches(valid); end
if (~isempty(reopatches)),reopatches = reopatches(reovalid); end
if (~isempty(ahspatches)),ahspatches = ahspatches(ahsvalid); end
if (~isempty(rhspatches)),rhspatches = rhspatches(rhsvalid); end
if (~isempty(prerunpatches)), prerunpatches = prerunpatches(prerunvalid); end
if (~isempty(postrunpatches)), postrunpatches = postrunpatches(postrunvalid); end

% 
% if (~exist('prerunpatches', 'var'))
%     prerunpatches = [];
% else
%     prerunpatches = prerunpatches(prerunvalid);
% end
% if (~exist('postrunpatches', 'var'))
%     postrunpatches = [];
% else
%     postrunpatches = postrunpatches(postrunvalid);
% end