classdef GlobalLookupTable < GlobalQuantity
    %UNTITLED5 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        highResDerivationMethod = @GlobalQuantity.oneDinterpolation_Clipped;
        
    end
    
    methods
        
       
        %{
        function data = eventTriggeredDataMatrix(glt, expt, trackpart, position, varargin)
            if (isa(expt, 'ExperimentSet'))
                data = [];
                for j = 1:length(expt.expt)
                    data = [data eventTriggeredDataMatrix(glt, expt.expt(j), trackpart, position, varargin{:})];
                end
                return
            end
            
            xField = glt.xField;
            varargin = assignApplicable(varargin);
        end
        %}
        
        function data = makeDataMatrixStretched(glt, starts, stops, npoints, varargin)
            %like makeDataMatrix, but instead each row running over
            %center+displacement axis, each row runs over
            %linspace(starts(j), stops(j), npoints)
            
           starts = starts(:);
           stops = stops(:);
           
           mult1 = repmat(linspace(0,1,npoints), size(starts));
           mult2 = 1-mult1;
           da = reshape(repmat(starts, [npoints 1]), size(mult1)).*mult2 + reshape(repmat(stops, [npoints 1]), size(mult2)).*mult1;
           data = reshape(interp1(glt.xData, glt.yData, da(:), 'linear', NaN), size(da));          
        end
        
        function [avg,datamatrix] = triggeredStretchedAverage (glt, starts, stops, npoints, varargin)
           datamatrix = makeDataMatrixStretched(glt, starts, stops, npoints, varargin{:});
           datamatrix2 = datamatrix;
           datamatrix2(~isfinite(datamatrix)) = 0; %replace unknown values with 0
           avg =  sum(datamatrix2)./sum(isfinite(datamatrix)); %average is sum of known values / number of known values
            
        end
        
        
        function data = makeDataMatrix(glt, centers, displacementAxis, varargin)
            % data = makeDataMatrix(glt, centers, displacementAxis, varargin) 
            % data is a matrix of values; each (row/column -- figure out
            % Natalie) represents interpolated values of the lookupTable at
            % various times. The times in (row/column) j are given by
            % centers(j) + displacementAxis
            %
            % e.g. if centers(3) is 104 and displacementAxis is -10:0.1:10,
            % then the 3rd (row/column) will represent the values of the
            % lookup table at times 94:0.1:114
            %
            % (behavior change) the value will be the average of the values
            % in the bin rather than the interpolated value at the bin
            % center
            %
            %glt < GlobalLookupTable
            %centers - M vector of center positions
            %displacmentAxis - N element list of offsets from center
            %position
            %varargin - not used currently
            
            centers = centers(:);
            displacementAxis = displacementAxis(:)';
            da = repmat(displacementAxis, size(centers));
            
            
            da = da + reshape(repmat(centers, size(displacementAxis')), size(da));
            
            
            dt = median(diff(displacementAxis));
            %08/12/2015 fixed bug where dt negative not treated correctly; 
            %could cause incorrect kernel fits, where (taxis is usually
            %descending)
            if (isfinite(dt) && abs(dt) > median(diff(glt.xData)))
                %     kk = ones([1 ceil(abs(dt)/median(diff(glt.xData)))]);
                %    kk = kk/(sum(kk));
                
                %adjusted 7/25/2016 to fix minor statistical aberration
                dtbins = abs(dt)/median(diff(glt.xData));
                kk = ones([1 ceil(dtbins)]);
                kk(1) = 1 - (ceil(dtbins)-dtbins)/2;
                kk(end) = kk(1);
                kk = kk/sum(kk);
                
                yd = conv(glt.yData, kk, 'valid');
                xc = interp1(1:length(glt.xData), glt.xData, (1:length(yd)) - 0.5 + length(kk)/2);
                
                %                 xe = (min(glt.xData)-dt/2):dt:(max(glt.xData)+dt/2);
                %                 [xc,yd] = meanyvsx(glt.xData, glt.yData, xe);
                data = reshape(interp1(xc, yd, da(:), 'linear', NaN), size(da));
            else
                if (islogical(glt.yData))
                    data = logical (reshape(interp1(glt.xData, double(glt.yData), da(:), 'nearest', 'extrap'), size(da)));
                else
                    data = reshape(interp1(glt.xData, glt.yData, da(:), 'linear', NaN), size(da));
                end
            end
            %data = reshape(glt.highResDerivationMethod(da, glt.xData, glt.yData), length(displacementAxis), [])';
        end
        
        function [avg,datamatrix] = triggeredAverage (glt, centers, displacementAxis, varargin)
            %[avg,datamatrix] = triggeredAverage (glt, centers, displacementAxis, varargin)
            % todo comment = triggered average just takes the average of
            % the data matrix, keeping track of out-of-range values
            
           datamatrix = makeDataMatrix(glt, centers, displacementAxis, varargin{:}); 
           datamatrix2 = datamatrix;
           datamatrix2(~isfinite(datamatrix)) = 0; %replace unknown values with 0
           avg =  sum(datamatrix2)./sum(isfinite(datamatrix)); %average is sum of known values / number of known values
            
        end
        
        function gq = toGlobalQuantity(glt, expt, addToExpt, varargin)
            %function gq = toGlobalQuantity(glt, expt, addToExpt, varargin)
            existsAndDefault('addToExpt', true);
            
            if (length(glt) > 1)
                for j = 1:length(glt)
                    gq(j) = toGlobalQuantity(glt(j),expt,addToExpt,varargin{:}); %#ok<AGROW>
                end
                return;
            end
            gq = GlobalQuantity();
            if (~strcmpi(glt.xField, 'eti'))
                gq = repmat(gq,0);
                return;
            end
            gq.xField = 'eti';
            gq.fieldname = glt.fieldname;
            gq.xData = unique(expt.elapsedTime(:)');
            gq.xData = gq.xData(isfinite(gq.xData));
            gq.yData = GlobalLookupTable.averageInPrecedingBin(gq.xData, glt.xData, glt.yData);
            if (addToExpt)
                expt.addGlobalQuantity(gq);
            end
        end
        function [u_exp, s_exp] = theoreticalSumSqStats (glt, tta_time, ntime_bins, nturns)
            % function [u_exp, s_exp] = theoreticalSumSqStats (glt, tta_time, ntime_bins, nturns)
            
            if (length(glt) > 1)
                if (length(nturns) == 1)
                    nturns = nturns * ones(size(glt));
                end
                u = zeros(size(glt));
                s = zeros(size(glt));
                for j = 1:length(glt)
                    [u(j),s(j)] = glt(j).theoreticalSumSqStats(tta_time, ntime_bins, nturns(j));
                end
                
                u_exp = sum(nturns.^2.*u)./sum(nturns)^2;
                s_exp = sqrt(sum(nturns.^3.*s.^2)./sum(nturns)^3);
                return;
            end
            
            vy = var(glt.yData);
            dt = median(diff(glt.xData));
            dtau = tta_time./(ntime_bins - 1);
            vx = vy .* dt./dtau;
            u_exp = vx.*ntime_bins./nturns;

            s_exp = sqrt(2*ntime_bins./nturns.^2.*vy^2.*(dt./dtau).^2);
        end
        function [u, s, vals] = shiftedSumSqStats (glt, centers, displacementAxis)
            
            trange = [min(centers) max(centers)];
            deltaT = max(displacementAxis) - min(displacementAxis);
            nshifts = floor(diff(trange)/deltaT);
            
            vals = zeros([1 nshifts]);
            
            for j = 1:nshifts
                newt = mod(centers + j*deltaT - trange(1), diff(trange)) + trange(1);
                vals(j) = sum(glt.triggeredAverage(newt, displacementAxis).^2);
            end
            
            u = mean(vals);
            s = std(vals);
        end
        
    end
    
    methods (Static)
         [ledLookupTable, ledGlobalQuantity] = createLedTableFromBitFile(expt, bitfilename, addToExpt)
         
         function yout = averageInPrecedingBin(xin, xData, yData) 
             %function yout = averageInPrecedingBin(xin, xData, yData)
             tp = false;
             if (size(xin,1) > size(xin, 2))
                 xin = xin';
                 tp = true;
             end
             [sortx,I] = sort(xin);[~,J] = sort(I);
             xbin = [2*sortx(1)-sortx(2) sortx]; %bug fixed 2/17 by MHG; don't know how it worked at all before???
             [~,yout] = meanyvsx(xData, double(yData), xbin);
             yout = yout(J);
             if (tp)
                 yout = yout';
             end
         end
             
         
    end
end

