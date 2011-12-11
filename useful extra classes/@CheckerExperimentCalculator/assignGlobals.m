function assignGlobals(cec, expt,varargin)
%reassigns all tracked global quantities to all tracks
%function assignGlobalQuantities(expt,varargin)
%expt.assignGlobalQuantities
%
%reassigns all tracked global quantities (see expt.globalQuantity)
%to all tracks
%'fields', {fieldnames} only assings globals with given field names

if (length(expt) > 1)
    for j = 1:length(expt)
        assignGlobals(cec, expt(j), varargin{:});
        disp ([num2str(j) '/' num2str(length(expt)) ' experiments assigned']);
    end
    return;
end

fieldnames = {};
varargin = assignApplicable(varargin);
ts0 = tic;

if (isempty(fieldnames))
    kinds = 1:length(cec.globalQuantities);
else
    [~,kinds] = intersect({cec.globalQuantities.fieldname}, fieldnames);
end
    

for k = kinds
    ts1 = tic;
    lasttime = toc(ts1); 
    gq = cec.globalQuantities(k);
    %{
    for j = 1:length(expt.track)
        gq.addQuantityToTrack(expt.track(j));
        if (toc(ts1) - lasttime > 60)
            disp ([num2str(j) ' - ' num2str(toc(ts1))]); 
            lasttime = toc(ts1);
        end
    end
    %}
    expt.addGlobalQuantity(gq);
    if (toc(ts0) > 60)
        disp ([num2str(k) ' / ' num2str(length(kinds)) ' assigned']);
        ts0 = tic;
    end
end
