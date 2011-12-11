classdef GlobalQuantity
    %GlobalQuantity specifies an additional field that is derived
    %from the position/time/etc. of a track and some additional information
    %that is the same for all tracks
    %
    %an example:  say the temperature varies with time, then you could
    %define a global quantity with
    %independent field 'eti' (interpolated time)
    %field name 'temperature'
    %xData (a sequence of times)
    %yData (a sequence of temperatures at those times)
    %derivation method: GlobalQuantity.oneDinterpolation
    
    properties
        xField = 'eti'; %name of the derived quantity to which the derivationMethod will be applied                                               
        fieldname = ''; %name of the field in which the global data will be stored
        xData = []; %xData, yData form pairs mapping the xField to the new field
        yData = []; %xData, yData form pairs mapping the xField to the new field
        derivationMethod = @GlobalQuantity.oneDinterpolation; %yfield = derivationMethod (xfield, xData, yData);
       
        %NOTES:
        %track.getDerivedQuantity(xField) must return valid data
        %derivationMethod is a function of the form
        %yfield = derivationMethod (xfield, xData, yData);
        %see oneDinterpolation and twoDinterpolation in this file for an example
     
    end
    
    methods
        function addQuantityToTrack (gq, track)
            %applies Global Quantity to Track using derivationMethod
            if (ischar (gq.xField))
                track.dq.(gq.fieldname) = gq.derivationMethod(track.getDerivedQuantity(gq.xField), gq.xData, gq.yData);
            else
                for j = 1:length(gq.xField)
                    xin.(gq.xField{j}) = track.getDerivedQuantity(gq.xField{j});                 
                end
                track.dq.(gq.fieldname) = gq.derivationMethod(xin, gq.xData, gq.yData);
            end
        end
        gqs = timeOnOffGQs (gq, ramptype, varargin);
    end
        
    
    methods (Static)
        function yout = oneDinterpolation(xin, xData, yData)
            % does a 1-D lookup table
            % function yout = oneDinterpolation(xin, xData, yData)

            %reshape y to have correct dimension
            sz = size(yData);
            n = length(xData);
            
            ind = find(sz == n);
            ind2 = find(sz ~= n);
            
            if (~any(ind))
                errmsg = ['yData and xData have incommensurate sizes: size(yData) = ' mat2str(sz) ...
                    ' and size xData = ' mat2str(size(xData))];
                error('GERSHOW:GQ01', errmsg);
            end
            
            order = [ind ind2];
            yData = permute(yData, order);
            size(yData);
            
            %make sure xin is a column vector
            if (size(xin,2) == 1)
                xin = xin';
            end
            s = warning('off','all');
            %interpolate
            yout = interp1 (double(xData), double(yData), xin, 'linear', NaN);
            warning(s);
            %reshape yout to have same shape as yData
         %   size(yout)
          %  yout = ipermute(yout, order);
           % size(yout)
        end
        function yout = oneDinterpolation_Clipped(xin, xData, yData)
            % does a 1-D lookup table with linear interpolation
            % function yout = oneDinterpolation(xin, xData, yData)
            % prior to interpolation, if xin < min(xData), xin is set to min(xData) and 
            % if xin > max(xData), xin is set to max(xData);
            
            %reshape y to have correct dimension
            sz = size(yData);
            n = length(xData);
            
            ind = find(sz == n);
            ind2 = find(sz ~= n);
            
            if (~any(ind))
                errmsg = ['yData and xData have incommensurate sizes: size(yData) = ' mat2str(sz) ...
                    ' and size xData = ' mat2str(size(xData))];
                error('GERSHOW:GQ01', errmsg);
            end
            
            order = [ind ind2];
            yData = permute(yData, order);
            size(yData);
            
            %make sure xin is a column vector
            if (size(xin,2) == 1)
                xin = xin';
            end
            s = warning('off','all');
            %interpolate
            xin(xin < min(xData)) = double(min(xData));
            xin(xin > max(xData)) = double(max(xData));
            yout = interp1 (double(xData), double(yData), xin, 'linear', NaN);
            warning(s);
            %reshape yout to have same shape as yData
         %   size(yout)
          %  yout = ipermute(yout, order);
           % size(yout)
        end
        
         function yout = twoDinterpolation(xin, xData, yData)
        % does a 2D lookup table
        % function yout = twoDinterpolation(xin, xData, yData)
        % xData is a NYxNXx2 matrix of x & y locations (x = xData(:,:,1)
        % and y = xData(:,:,2))
        % or xData is a struct with fields x and y
        % yData is a (NYxNXxk) matrix of values
        % xin is a 2XM vector of x&y locations
        % yout is a kXM vector
        
            if (isstruct(xData))
                x = xData.x;
                y = xData.y;
            else
                if (size(xData,3) > 1)
                    x = xData(:,:,1);
                    y = xData(:,:,2);
                else
                    x = xData(1,:);
                    y = xData(2,:);
                end
            end
            x = double(x);
            y = double(y);
            isl = islogical(yData);
            yData = double(yData);
            yout = NaN([size(yData,3) size(xin,2)]);
            inds = isfinite(xin(1,:)) & isfinite(xin(2,:));
            
            for k = 1:size(yData,3)
                yout(k,inds) = interp2(x,y,yData(:,:,k), xin(1,inds), xin(2,inds), 'linear');
            end
            if (isl)
                yout = logical(round(yout));
            end
         end
         function yout = interpAngleIm(xin, xdata, ydata)
             %interpolates 2D angle field correctly
             u = GlobalQuantity.twoDinterpolation(xin, xdata, cos(ydata));
             v = GlobalQuantity.twoDinterpolation(xin, xdata, sin(ydata));
             yout = atan2(v,u);
         end
         yout = tri2Dinterpolation(xin, xData, yData);
         gasval = timeVaryingGasDerivation (xin, xdata, ydata);
         gasval = timeVaryingGasDerivationNearest (xin, xdata, ydata);
         yout = interpLogicalIm(xin, xdata, ydata);
         yout = oneDinterpolationNearest(xin, xData, yData)
    end
    
end

