function [qv, datamatrix] = averageFromSubField(expt, subfield, field, centerpos, offsetinds, varargin)
%TODO: update comments
%gathers fieldname from sections of track that are part of subfield
%function qvec = gatherFromSubField(expt, subfield, fieldname, varargin)
%
%gets derived quantity fieldname from subfield by calling
%track.getSubFieldDQ(subfield, fieldname, varargin)
%
%see Track.getSubFieldDQ and 
%    TrackPart.getDerivedQuantity for optional params
%later we should update to use full syntax of track.getSubFieldDQ, but for
%now let's just get working

if (length(expt) > 1)
    [~,datamatrix] =  averageFromSubField(expt(1), subfield, field, centerpos, offsetinds, varargin{:});
    for j = 2:length(expt)
        [~,dm] = averageFromSubField(expt(j), subfield, field, centerpos, offsetinds, varargin{:});
        datamatrix = [datamatrix;dm]; %#ok<AGROW>
    end
    qv = mean(datamatrix);
    return;
end
            
[qv, datamatrix] = averageDerivedQuantity([expt.track.(subfield)], field, centerpos, offsetinds, varargin{:});
 