function qvec = gatherFromSubField(expt, subfield, fieldname, varargin)
%gathers fieldname from sections of track that are part of subfield
%function qvec = gatherFromSubField(expt, subfield, fieldname, varargin)
%
%gets derived quantity fieldname from subfield by calling
%track.getSubFieldDQ(subfield, fieldname, varargin)
%
%see Track.getSubFieldDQ and 
%    TrackPart.getDerivedQuantity for optional params

 qvec = [];
 for j = 1:length(expt.track)
     qvec = [qvec expt.track(j).getSubFieldDQ(subfield, fieldname, false, varargin{:})];
     %{
     if (isempty(qvec))
         break;
     end
     %} 
     %removed by MHG 8/23/2010 - don't know what purpose this served
 end
 
 