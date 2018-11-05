classdef BehaviorTriggeredData
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        fname;
        es; %experiment statistics
        gq;
        glt;
        hs;
        reo;
        run;
        turn;
        acchs;
        rejhs;
        pause;
        all;
    end
    methods(Static)
       btdstruct = analyzeBTDDirectory_onesource (btdstruct, gqname, varargin);
       btdstruct = loadBTDDirectory (basedir, varargin);
       btdstruct = analyzeBTDDirectory_twosource (btdstruct, stim1name, stim2name, varargin);
       btdstruct = analyzeVarianceSwitchingExperiment(btdstruct, opstruct);
       btdstruct = analyzeVarianceSwitchingExperiment_Gepner(btdstruct, opstruct, qfactor);
       %problemDescription = createProblemStructForRateFunFitting (var, timeField, polynomialDegree, period, trange);
       problemDescription = createProblemStructForRateFunFitting (vs, timeField, polynomialDegree, trange, exprange, deltaT);
       problemDescription = simulateVarSwitchFromPrevFit (btdstruct, osfit, oldpd, alphaLow, alphaHigh, nlarvae);
       btdstruct = prepVarianceSwitchingAnalysis(btdstruct, opstruct);
       btdstruct = prepVarianceSwitchingAnalysis_Gepner(btdstruct, opstruct);
       btdstruct = prepVarianceSwitchingAnalysis_Wolk(btdstruct, opstruct);
       ah = summaryFigure(bdstruct, varargin);    
       
       function [dkl12, tcent, bit_rate, tstart] = klDivergenceBetweenDM (taxis, dm1, dm2, tbinEdges, varargin)
       %function [dkl12, tcent] = klDivergenceBetweenDM (taxis, dm1, dm2, tbinEdges, varargin)
       %function [dkl12, tcent] = klDivergenceBetweenDM (taxis, dm1, dm2, binWidth, varargin)
       %dkl = dkl(x1||x2)
           method = 'gauss';
           varargin = assignApplicable(varargin);
           
           if (size(tbinEdges,1) == 1)
               if (size(tbinEdges,2) == 1) %it's a bin width rather than edges
                    bw = tbinEdges;
                    tbinEdges = taxis(taxis < max(taxis - bw));
                    tbinEdges(2,:) = tbinEdges(1,:) + bw;
               else
                   tbinEdges(2,:) = [tbinEdges(2:end) tbinEdges(end)+eps];
               end
           end
           
           %[~,bin] = histc(taxis, tbinEdges);
           dkl12 = zeros(1, size(tbinEdges,2)-1);
           tcent = dkl12;
           for j = 1:(size(tbinEdges,2)-1)
               inds = taxis >= tbinEdges(1,j) & taxis < tbinEdges(2,j);
               tcent(j) = mean(taxis(inds));
               x1 = dm1(:, inds);
               x2 = dm2(:, inds);
               x1 = x1(all(isfinite(x1),2), :); %#ok<CPROP> all is also a property; function meant here
               x2 = x2(all(isfinite(x2),2), :); %#ok<CPROP> all is also a property; function meant here
               if (isempty(x1) || isempty(x2))
                   continue;
               end
               switch (lower(method))
                   case 'gauss'
                       dkl12(j) = DKL_gauss_est(x1, x2);
                   case 'knnti'
                       if (j == 1)
                           co = DKL_kNN_kiTi_initialization(1, {'kNNmethod', 'knnFP1'});
                       end
                       dkl12(j) = DKL_kNN_kiTi_estimation(x1', x2', co);
                   otherwise
                       error(['unknown method: ' method]);
               end
           end
           bit_rate = (dkl12/log(2)) ./ diff(tbinEdges(:,1:end-1));
           tstart = tbinEdges(1,1:end-1);
       end
       
       
       
       function [mi, tcent, bit_rate, tstart] = miDiscreteVsTime (taxis, dms, tbinEdges, varargin)
       %function [mi, tcent] = miDiscreteVsTime (taxis, dms, tbinEdges, varargin)
       %dms should be a cell of datamatrices, each corresponding to a
       %discrete state
            dmall = vertcat(dms{:});
            ntotal = size(dmall, 1);
            [mi,tcent, bit_rate, tstart] = BehaviorTriggeredData.klDivergenceBetweenDM(taxis, dms{1}, dmall, tbinEdges, varargin{:});
            mi = mi*size(dms{1},1)/ntotal;
            bit_rate = bit_rate*size(dms{1},1)/ntotal;
            for j = 2:length(dms)
                [mitemp, ~,brtemp] = BehaviorTriggeredData.klDivergenceBetweenDM(taxis, dms{j}, dmall, tbinEdges, varargin{:});
                mi = mi + size(dms{j},1)/ntotal*mitemp;
                bit_rate = bit_rate + size(dms{j},1)/ntotal*brtemp;               
            end
       
       end
       
       
       
    end
    
    methods
        function btd = BehaviorTriggeredData(expt, varargin)

                %note -- this is fine syntax, but for some reason causes MATLAB
                %to run exceedingly slow.  run loop elsewhere
    %         if (length(expt) > 1)
    %             ts1 = tic;
    %             for j = 1:length(expt)
    %                 btd(j) = BehaviorTriggeredData(expt(j), varargin{:}); %#ok<AGROW>
    %             end
    %             disp ([num2str(j) '/' num2str(length(expt)) ' ' num2str(toc(ts1),3)]); 
    %             return;
    %         end

            %Gather all the fields from the experiment
            fields = {};
            
            varargin = assignApplicable(varargin);

            fields = union(fields, {'eti', 'speed', 'spineTheta', 'curv', 'deltatheta', 'vel_dp', 'lrdtheta', 'vhead', 'vheadperp', 'acc'});


            btd.fname = expt.fname;
            btd.glt = expt.globalLookupTable;
            %btd.gq = expt.globalQuantity;
            for j = 1:length(btd.glt)
                fields = union(fields, btd.glt(j).xField);
            end
            
            btd.es = expt.getExperimentStatistics();
% 
%             tx = min(expt.gatherField('eti')):(100*expt.dr.interptime):max(expt.gatherField('eti'));
%             btd.numAnimals = max(

            btd.hs.tnum = expt.gatherFromSubField('headSwing', 'trackNum');
            btd.reo.tnum = expt.gatherFromSubField('reorientation', 'trackNum');
            btd.run.tnum = expt.gatherFromSubField('run', 'trackNum');

            rs = expt.gatherSubField('run', 'startTheta');
            re = expt.gatherSubField('run', 'endTheta');
            rt = expt.gatherSubField('run', 'runTime');
            rmt = expt.gatherSubField('run', 'meanTheta');

            dt = diff(unwrap([rs;re]));

            btd.run.startTheta = rs;
            btd.run.endTheta = re;
            btd.run.time = rt;
            btd.run.meanTheta = rmt;
            btd.run.deltaTheta = dt;
            btd.run.length = expt.gatherSubField('run', 'pathLength');
            btd.run.euclidLength = expt.gatherSubField('run', 'euclidLength');

            btd.reo.numhs = expt.gatherSubField('reorientation', 'numHS');
            btd.reo.nextdir = expt.gatherSubField('reorientation', 'nextDir');
            btd.reo.prevdir = expt.gatherSubField('reorientation', 'prevDir');
            btd.reo.dtheta = diff(unwrap([btd.reo.prevdir;btd.reo.nextdir]));
            
%             r = expt.gatherField('reorientation');
%             dr = [expt.track.dr]; 
%             preReoTime = 5; %seconds
%             npts = ceil(preReoTime/median([dr.interpTime]));
%             btd.reo.periPhaseMatrix = zeros(length(r), npts);
%             btd.reo.periPhaseMatrixAdj = zeros(length(r), npts);
%             btd.reo.timeMatrix = zeros(length(r), npts);
%             px = 0:(-1):(1-npts);
%             for j = 1:length(r)
%                 pp = r(j).getDerivedQuantity('periPhase', 'position', 'startInd', 'posoffset', px);
%                 dpp = diff(pp);
%                 ind1 = find(dpp < 0, 1, 'first');
%                 ind2 = find(dpp(ind1:end) > 4, 1, 'first') + ind1;
%                 ppadj = unwrap(pp);
%                 ppadj = ppadj - ppadj(ind2);
%                 btd.reo.periPhaseMatrix(j,:) = pp;
%                 btd.reo.periPhaseMatrixAdj(j,:) = ppadj;
%                 btd.reo.timeMatrix(j,:) = r(j).getDerivedQuantity('eti', 'position', 'startInd', 'posoffset', px);
%             end
%                 
            
            btd.hs.taildir = expt.gatherSubField('headSwing', 'tailDir');
            btd.hs.headdir = expt.gatherSubField('headSwing', 'headDir');
            btd.hs.accepted = expt.gatherSubField('headSwing', 'accepted');
            btd.hs.sign = expt.gatherSubField('headSwing', 'sign');
            btd.hs.htv = logical(expt.gatherSubField('headSwing', 'valid'));
            btd.hs.maxtheta = expt.gatherSubField('headSwing', 'maxTheta');
            btd.hs.hsnum = expt.gatherSubField('headSwing', 'num');

            for j = 1:length(fields)
                btd.run.(['start_' fields{j}]) = expt.gatherFromSubField('run', fields{j}, 'position', 'start');
                btd.run.(['end_' fields{j}]) = expt.gatherFromSubField('run', fields{j}, 'position', 'end');
                btd.reo.(['start_' fields{j}]) = expt.gatherFromSubField('reorientation', fields{j}, 'position', 'start');
                btd.reo.(['end_' fields{j}]) = expt.gatherFromSubField('reorientation', fields{j}, 'position', 'end');
                btd.hs.(['start_' fields{j}]) = expt.gatherFromSubField('headSwing', fields{j}, 'position', 'start');
                btd.hs.(['end_' fields{j}]) = expt.gatherFromSubField('headSwing', fields{j}, 'position', 'end');
                btd.hs.(['atmax_' fields{j}]) = expt.gatherFromSubField('headSwing', fields{j}, 'position', 'atmax');
                btd.all.(fields{j}) = expt.gatherField(fields{j});
            end
            btd.all.tnum = expt.gatherField('trackNum', 'expandToInterped', true);
            btd.all.isrun = expt.gatherField('isrun');

            %Separate the turns and pauses from the reorientations
            inds = [btd.reo.numhs] > 0;
            fn = fieldnames(btd.reo);
            for j = 1:length(fn)
                btd.turn.(fn{j}) = btd.reo.(fn{j})(:,inds);
            end

            for j = 1:length(fn)
                btd.pause.(fn{j}) = btd.reo.(fn{j})(:,~inds);
            end

            %Separate the accepted and rejected head sweeps
            inds = btd.hs.accepted;
            fn = fieldnames(btd.hs);
            for j = 1:length(fn)
                btd.acchs.(fn{j}) = btd.hs.(fn{j})(inds);
            end
            for j = 1:length(fn)
                btd.rejhs.(fn{j}) = btd.hs.(fn{j})(~inds);
            end
            
        end%function
        
        function btd = addTonToffGQs (btd, fieldname, ramptype, varargin)
            existsAndDefault('ramptype', 'square');
            for j = 1:length(btd)
                ind = btd(j).findField(fieldname);
                if (ind <= 0)
                    error ([fieldname ' not a global lookup table in btd']);
                end
                
                btd(j) = btd(j).addGlobalLookupTable(btd(j).glt(ind).timeOnOffGQs(ramptype, varargin{:}));
                
            end
        end
        
        function btd = addVarianceGQs (btd, fieldname, ramptype, varargin)
            existsAndDefault('ramptype', 'square');
            for j = 1:length(btd)
%                 if (btd(j).findField([fieldname '_var_ton']) > 0)
%                     continue;
%                 end
                
                ind = btd(j).findField(fieldname);
                if (ind <= 0)
                    error ([fieldname ' not a global lookup table in btd']);
                end
                
                btd(j) = btd(j).addGlobalLookupTable(btd(j).glt(ind).timeOnOffGQsVariance(ramptype, varargin{:}));
                
            end
        end

        function centers = getBTDMCenters (btd, tp, position, gltname, varargin) 
            gltind = btd.findField(gltname);
            conditions = [];
            varargin = assignApplicable(varargin);
            
            if (strcmpi(tp, 'all'))
                centers = btd.(tp).(btd.glt(gltind).xField);
            else
                centers = btd.(tp).([position '_' btd.glt(gltind).xField]);
            end
            if (~isempty(conditions))
                if (~isfield(conditions, 'name') || ~isfield(conditions, 'validop'))
                    warning ('conditions must have field name and validop; valid = validop(btd.tp.name)');
                else
                    valid = true(size(centers));
                    for j = 1:length(conditions)
                        if (isfield(btd.(tp), conditions(j).name))
                            valid = valid & conditions(j).validop(btd.(tp).(conditions(j).name));
                        else
                            vind = btd.findField(conditions(j).name);
                            if (vind > 0)
                                
                                if (strcmpi(tp, 'all'))
                                    c2 = btd.(tp).(btd.glt(vind).xField);
                                else
                                    c2 = btd.(tp).([position '_' btd.glt(vind).xField]);
                                end
                                valid = valid & conditions(j).validop(btd.glt(vind).highResDerivationMethod(c2, btd.glt(vind).xData, btd.glt(vind).yData));
                            else
                                warning ([conditions(j).name ' does not name a field in ' tp ' or a glt field']);
                            end
                        end
                    end
                    centers = centers(valid);
                end
            end
        end
        
        %Takes data in btd for a given track part and reformats it as a dataMatrix
        function [datamatrix, stimnum, btdnum] = behaviorTriggeredDataMatrix (btd, tp, position, gltname, displacementAxis, varargin)
            %function [datamatrix, stimnum] = behaviorTriggeredDataMatrix (btd, tp, position, gltname, displacementAxis, varargin)
            %tp=track part (run,reo,etc)
            %position=position within the track part (start, end, max,
            %etc.)
            
            if (length(btd) > 1)
                [datamatrix, stimnum, btdnum] =  behaviorTriggeredDataMatrix (btd(1), tp, position, gltname, displacementAxis, varargin{:});
                for j = 2:length(btd)
                    if (isempty(btd(j).glt))
                        continue;
                    end
                    [dm,sn,bn] =  behaviorTriggeredDataMatrix(btd(j), tp, position, gltname, displacementAxis, varargin{:});
                 	datamatrix = [datamatrix;dm]; %#ok<AGROW>
                    stimnum = [stimnum;sn];%#ok<AGROW>
                    btdnum = [btdnum;j*ones(size(bn))]; %#ok<AGROW> %btdnum is just 1 number per line
                end
                return;
            end
            
            if (iscell(gltname))
                [datamatrix,~,btdnum] = behaviorTriggeredDataMatrix(btd, tp, position, gltname{1}, displacementAxis, varargin{:});
                stimnum = ones(size(datamatrix));
                for k = 2:length(gltname)
                    dm = behaviorTriggeredDataMatrix(btd, tp, position, gltname{k} , displacementAxis, varargin{:});
                    datamatrix = [datamatrix dm]; %#ok<AGROW>
                    stimnum = [stimnum k*ones(size(dm))]; %#ok<AGROW>
                    %{
                    datamatrix = [datamatrix behaviorTriggeredDataMatrix(btd, tp, position, gltname{k} , displacementAxis, varargin{:})]; %#ok<AGROW>
                    stimnum = [stimnum k*ones(size(datamatrix))]; %#ok<AGROW>
                    %}
                end
                return;
            end
            btdnum = [];
            stimnum = [];
            conditions = [];
            shift = [];
            %conditions.name = 'hsnum';
            %conditions.validop = @(x) x == 1;
            %this will get first headsweep
           
            varargin = assignApplicable(varargin);
            
            
            %Make a datamatrix from the glt
            if (isempty(btd.glt))
                datamatrix = [];
                return;
            end
            try
                gltind = find(strcmpi(gltname, {btd.glt.fieldname}),1,'first');
            catch
                datamatrix = [];
                return;
            end
            if (isempty(gltind))
                datamatrix = [];
                return;
            end
            try
                centers = btd.getBTDMCenters(tp, position, gltname, 'conditions', conditions, varargin{:});
                if (~isempty(shift))
                    trange = [min(centers) max(centers)];
                    centers = mod(centers + shift - trange(1), diff(trange)) + trange(1);
                end
                datamatrix = btd.glt(gltind).makeDataMatrix(centers, displacementAxis, varargin{:});
%                 if (strcmpi(tp, 'all'))
%                     centers = btd.(tp).(btd.glt(gltind).xField);
%                 else
%                     centers = btd.(tp).([position '_' btd.glt(gltind).xField]);
%                 end
%                 if (~isempty(conditions))
%                     if (~isfield(conditions, 'name') || ~isfield(conditions, 'validop'))
%                         warning ('conditions must have field name and validop; valid = validop(btd.tp.name)');
%                     else
%                         valid = true(size(centers));
%                         for j = 1:length(conditions)
%                             if (isfield(btd.(tp), conditions(j).name))                             
%                                 valid = valid & conditions(j).validop(btd.(tp).(conditions(j).name));
%                             else
%                                 vind = btd.findField(conditions(j).name);
%                                 if (vind > 0)
%                                     
%                                     if (strcmpi(tp, 'all'))
%                                         c2 = btd.(tp).(btd.glt(vind).xField);
%                                     else
%                                         c2 = btd.(tp).([position '_' btd.glt(vind).xField]);
%                                     end
%                                     valid = valid & conditions(j).validop(btd.glt(vind).highResDerivationMethod(c2, btd.glt(vind).xData, btd.glt(vind).yData));
%                                 else
%                                     warning ([conditions(j).name ' does not name a field in ' tp ' or a glt field']);
%                                 end
%                             end
%                         end
%                         centers = centers(valid);
%                     end
%                 end
               
                
            catch me
                disp(me.getReport());
                datamatrix = [];
            end
            stimnum = ones(size(datamatrix));
            btdnum = ones([size(datamatrix,1) 1]);
        end
        
        %Takes data in btd for a given track part and reformats it as a dataMatrix
        function [datamatrix, stimnum] = behaviorTriggeredDataMatrixCompensator (btd, tp, position, gltname, gltnameCompensation, displacementAxis, varargin)
            %function [datamatrix, stimnum] = behaviorTriggeredDataMatrixCompensator (btd, tp, position, gltname, displacementAxis, varargin)
            %tp=track part (run,reo,etc)
            %
            %finds value of gltnameCompensation at tp.position, subject to conditions then
            %randomly picks nreps (default = 10) other time points with same value of
            %gltnameCompensation (+/- max(epsilon_val, epsilon_pct*value), default epsilon_pct = 0.01, epsilon_val = 1)
            %returns data matrix for these centers
            %
            %additional args: 'epsilon_pct', 0.01, 'epsilon_val', 1; 'nreps', 10
            
            if (length(btd) > 1)
                [datamatrix, stimnum] =  behaviorTriggeredDataMatrixCompensator (btd(1), tp, position, gltname, gltnameCompensation, displacementAxis, varargin{:});
                for j = 2:length(btd)
                    if (isempty(btd(j).glt))
                        continue;
                    end
                    [dm,sn] =  behaviorTriggeredDataMatrixCompensator(btd(j), tp, position, gltname, gltnameCompensation, displacementAxis, varargin{:});
                 	datamatrix = [datamatrix;dm]; %#ok<AGROW>
                    stimnum = [stimnum;sn];%#ok<AGROW>
                end
                return;
            end
            
            if (iscell(gltname))
                datamatrix = behaviorTriggeredDataMatrixCompensator(btd, tp, position, gltname{1}, gltnameCompensation, displacementAxis, varargin{:});
                stimnum = ones(size(datamatrix));
                for k = 2:length(gltname)
                    dm = behaviorTriggeredDataMatrixCompensator(btd, tp, position, gltname{k}, gltnameCompensation, displacementAxis, varargin{:});
                    datamatrix = [datamatrix dm]; %#ok<AGROW>
                    stimnum = [stimnum k*ones(size(dm))]; %#ok<AGROW>
                    %{
                    datamatrix = [datamatrix behaviorTriggeredDataMatrix(btd, tp, position, gltname{k} , displacementAxis, varargin{:})]; %#ok<AGROW>
                    stimnum = [stimnum k*ones(size(datamatrix))]; %#ok<AGROW>
                    %}
                end
                return;
            end
            nreps = 10;
            epsilon_pct = 0.01;
            epsilon_val = 1;
            varargin = assignApplicable(varargin);
            if (isempty(btd.glt))
                datamatrix = [];
                stimnum = [];
                return;
            end
            compval = behaviorTriggeredDataMatrix(btd, tp, position, gltnameCompensation, 0, varargin{:});
            compval = compval(isfinite(compval));
            if (isempty(compval))
                datamatrix = [];
                stimnum = [];
                return;
            end
            compind = btd.findField(gltnameCompensation);
            gltind = btd.findField(gltname);
            if (compind <= 0 || gltind <= 0)
                datamatrix = [];
                stimnum = [];
                return;
            end
            
            
            %Make a datamatrix from the glt
            xdata = btd.glt(compind).xData;
            ydata = btd.glt(compind).yData;
            centers = zeros(length(compval), nreps);
            for j = 1:length(compval)
               inds = find(abs(compval(j)-ydata) <= max(epsilon_val, epsilon_pct*compval(j)));
               centers(j,:) = xdata(inds(ceil(rand([1 nreps])*length(inds))));  
            end
            
            datamatrix = btd.glt(gltind).makeDataMatrix(centers(:), displacementAxis, varargin{:});
            stimnum = ones(size(datamatrix));
        end
        
        function ind = findField(btd, fieldname)
            %function ind = findField(btd, fieldname)
            %ind = 0 on failure, otherwise 
            %strcmpi(btd.glt(ind).fieldname, fieldname) is true
             if (length(btd) > 1)
                 ind = zeros(size(btd));
                 for j = 1:length(btd)
                     ind(j) = btd(j).findField(fieldname);
                 end
                 return;
             end
             if (isempty(btd) || isempty(btd.glt)) %|| ~isa(btd.glt, 'GlobalQuantity'))
                 ind = 0;
                 return;
             end
             fn = {btd.glt.fieldname};
             ind = find(strcmpi(fieldname, fn), 1, 'first');
             if (isempty(ind))
                 ind = 0;
                 %warning ('BTD:findField', [fieldname ' not a valid glt name']);
             end
             
        end
        function glts = gatherGlt (btd, fieldname)
            inds = btd.findField(fieldname);
            if (any (inds <= 0))
                glts = [];
                return;
            end
            for j = 1:length(inds)
                glts(j) = btd(j).glt(inds(j)); %#ok<AGROW>
            end
        end
        
        
        function btd = addGlobalLookupTable (btd, glt)
            %function btd = addGlobalLookupTable (btd, glt)
            for j = 1:length(glt)
                ind = btd.findField(glt(j).fieldname);
                if (ind <= 0)
                    btd.glt = [btd.glt glt(j)];
                else
                    btd.glt(ind) = glt(j);
                end
            end
        end
        
        %Takes data in btd for a given track part and reformats it as a stretched dataMatrix
        function datamatrix = behaviorTriggeredDataMatrixStretched (btd, tp, positionStart, positionEnd, gltname, npoints, varargin)
            %function datamatrix = behaviorTriggeredDataMatrix (btd, tp, position, gltname, displacementAxis, varargin)
            disp ('warning: not sure this works correctly -- should be checked more carefully');
            
            if (length(btd) > 1)
                datamatrix =  behaviorTriggeredDataMatrixStretched (btd(1), tp, positionStart, positionEnd, gltname, npoints, varargin{:});
                for j = 2:length(btd)
                    dm =  behaviorTriggeredDataMatrixStretched(btd(j), tp, positionStart, positionEnd, gltname, npoints, varargin{:});
                 	datamatrix = [datamatrix;dm]; %#ok<AGROW>
                end
                return;
            end
            
            conditions = [];
            %conditions.name = 'hsnum';
            %conditions.validop = @(x) x == 1;
            %this will get first headsweep
            varargin = assignApplicable(varargin);
            
            if (isempty(btd.glt))
                datamatrix = [];
                return;
            end
            try
                gltind = find(strcmpi(gltname, {btd.glt.fieldname}),1,'first');
            catch
                datamatrix = [];
                return;
            end
            if (isempty(gltind))
                datamatrix = [];
                return;
            end
           
            try
                starts = btd.(tp).([positionStart '_' btd.glt(gltind).xField]);
                stops = btd.(tp).([positionEnd '_' btd.glt(gltind).xField]);
                if (~isempty(conditions))
                    if (~isfield(conditions, 'name') || ~isfield(conditions, 'validop'))
                        warning ('conditions must have field name and validop; valid = validop(btd.tp.name)');
                    else
                        valid = true(size(starts));
                        for j = 1:length(conditions)
                            valid = valid & conditions(j).validop(btd.(tp).(conditions(j).name));
                        end
                        starts = starts(valid);
                        stops = stops(valid);
                    end
                end
                
                datamatrix = btd.glt(gltind).makeDataMatrixStretched(starts, stops, npoints, varargin{:});
            catch
                datamatrix = [];
            end
        end
        
 
        %Takes the average of the behavior triggered data, returning the
        %average and the data matrix
        function [bta, dm, stimnum, sem, sd] = behaviorTriggeredAverage (btd, tp, position, gltname, displacementAxis, varargin)
            % function [bta, dm] = behaviorTriggeredAverage (btd, tp, position, gltname, displacementAxis, varargin)
            [dm,stimnum] = btd.behaviorTriggeredDataMatrix(tp, position, gltname, displacementAxis, varargin{:});
            dm2 = dm; dm2(~isfinite(dm2)) = 0;
            bta = sum(dm2)./sum(isfinite(dm));
            
            
            stimnum = mean(stimnum);
            
            if (nargout > 3)
                btmat = repmat(bta, [size(dm,1), 1]);
                
                dm2 = dm - btmat; dm2(~isfinite(dm2)) = 0;
                sd = sqrt(sum(dm2.^2)./sum(isfinite(dm)));
                sem = sd./sqrt(sum(isfinite(dm)));
            end
            
        end
        %Takes the average of the behavior triggered data, returning the
        %average and the data matrix
        function [bta, dm, stimnum, sem, sd] = behaviorTriggeredAverageCompensator (btd, tp, position, gltname, gltnameCompensation, displacementAxis, varargin)
            % function [bta, dm] = behaviorTriggeredAverage (btd, tp, position, gltname, displacementAxis, varargin)
            [dm,stimnum] = btd.behaviorTriggeredDataMatrixCompensator(tp, position, gltname, gltnameCompensation, displacementAxis, varargin{:});
            dm2 = dm; dm2(~isfinite(dm2)) = 0;
            bta = sum(dm2)./sum(isfinite(dm));
            
            
            stimnum = mean(stimnum);
            
            if (nargout > 3)
                btmat = repmat(bta, [size(dm,1), 1]);
                dm2 = dm - btmat; dm2(~isfinite(dm2)) = 0;
                sd = sqrt(sum(dm2.^2)./sum(isfinite(dm)));
                sem = sd./sqrt(sum(isfinite(dm)));
            end
            
        end
        function [convkernel,btd] = createBTAKernel (btd, tp, position, gltname, t0, dt, varargin)
            %function convkernel = createBTAKernel (btd, tp, position, gltname, displacementAxis, varargin)
            %note: t0 should be less than 0; if t0 > 0, we assume you mean t0 < 0 and invert
            %kernel is created at points t0:dt:0
            %
            %kernel is normalized so that sum(convkernel.^2) = 1;

            newFieldName = {};
            abbott = false;
            loadKernel = false;
            varargin = assignApplicable(varargin);
            
            displacementAxis = -(0:dt:abs(t0));

            [bta,~,stimnum] = behaviorTriggeredAverage (btd, tp, position, gltname, displacementAxis, varargin{:});
            convkernel = zeros(max(stimnum),length(displacementAxis));
          %  alldata = btd.behaviorTriggeredDataMatrix('all', [], gltname,0);
            for j = 1:max(stimnum)
                if (abbott)
                    t = displacementAxis;
                    if (loadKernel)
                        load('E:\Variance Adaptation\MID_pillcbf1_ber17_51_120_diff.mat');  %load kernel from file.  loaded kernel must match displacementAxis -(0:dt:abs(t0))
                        btaj = MID_pillcbf1_ber17_51_120_diff';    %load kernel from file.  loaded kernel must match displacementAxis -(0:dt:abs(t0))
                    else
                        btaj = bta(stimnum == j);
                    end

                    abbottfun = @(x,t) (x(3)*exp(-x(2)*t)/(x(1)-x(2)) - x(5)*exp(-x(4)*t)/(x(1)-x(4)) + (x(3)*(x(4)-x(1)) - x(5)*(x(2)-x(1)))*exp(-x(1)*t)/((x(1)-x(2))*(x(1)-x(4)))) .* (t <= 0);
                    
                    
                    [~,I] = max(abs(btaj));
                    mv = btaj(I);
                    x0 = [-1.5 -1 mv -0.5 mv*1.5];
                    
                    op = optimset('lsqcurvefit');
                    op.MaxFunEvals = 5E4;
                    op.MaxIter = 5000;
                    
                   % op.Display = 'none';
                    
                    x = lsqcurvefit(abbottfun, x0, displacementAxis, btaj,[],[],op);
                    convkernel(j,:) = abbottfun(x,displacementAxis);
                else
                    if (loadKernel)
                        load('E:\Variance Adaptation\MID_pillcbf1_ber17_51_120_diff.mat');
                        ck = MID_pillcbf1_ber17_51_120_diff';
                    else
                        ck = bta(stimnum == j);
                    end
                                        
                    ck = csaps(1:length(ck), ck, dt, 1:length(ck));
                    convkernel(j,:) = ck./sqrt(sum(ck.^2));%/std(alldata(isfinite(alldata(:,j)),j));
                end
            end
            if (~iscell(gltname))
                gltname = {gltname};
            end
            if (~iscell(newFieldName))
                if (iscell(gltname) && length(gltname) > 1)
                    return;
                end
                newFieldName = {newFieldName};
                
            end
            for j = 1:length(newFieldName)
                if (ischar(newFieldName{j}) && ~isempty(newFieldName{j}))
                    btd = btd.addConvolvedFields (gltname{j}, newFieldName{j}, convkernel(j,:), dt, 'scaleToSqr', true, varargin{:});
                end
            end
            end
            
        function [stim, sp, delta_t] = getPillowData (btd, tp, position, gltname, delta_t, varargin) 
            %function [stim, sp] = getPillowData (btd, tp, position, gltname) 
            %gets data for use in pillow routine simpleSTC
            
            conditions = [];
            existsAndDefault('delta_t', []);
            
            if (length(btd) > 1)
                [stim, sp, delta_t] = getPillowData(btd(1), tp, position, gltname, delta_t, varargin{:});
                for j = 2:length(btd)
                    [s,p, delta_t] = getPillowData(btd(j), tp, position, gltname, delta_t, varargin{:});
                    p = p + length(stim);
                    stim = [stim;s]; %#ok<AGROW>
                    sp = [sp p]; %#ok<AGROW>
                end
                return;
            end
            xt = [];
            xte = [];
            if (~iscell(gltname))
                gltname = {gltname};
            end
            for k = 1:length(gltname)
                ind = btd.findField(gltname{k});
                if (ind <= 0)
                    warning ('btd:getPillowData', [gltname{k} ' not a valid field']);
                    stim = [];
                    sp = [];
                    return;
                end
                lt = btd.glt(ind);
                if (~any(strcmpi(lt.xField, {'eti', 'et'}))) % add other timing fields here if needed
                    warning ('btd:getPillowData', [gltname{k} ' xField may not be time']);
                end

                if (isempty(delta_t))
                    delta_t = median(diff(lt.xData));
                end
                if (isempty(xt)) 
                    xt = min(lt.xData):(delta_t):max(lt.xData);
                    xte = min(lt.xData-delta_t/2):delta_t:max(lt.xData+delta_t/2);
                end
    %            stim = interp1(lt.xData, lt.yData, xt, 'linear');
                [~,stim(:,k)] = meanyvsx(lt.xData, lt.yData, xte);
                %whos stim
            end
            if (strcmpi(tp, 'all'))
                centers = btd.(tp).(lt.xField);
            else
                centers = btd.(tp).([position '_' lt.xField]);
            end
            
             if (~isempty(conditions))
                    if (~isfield(conditions, 'name') || ~isfield(conditions, 'validop'))
                        warning ('conditions must have field name and validop; valid = validop(btd.tp.name)');
                    else
                        valid = true(size(centers));
                        for j = 1:length(conditions)
                            if (isfield(btd.(tp), conditions(j).name))                             
                                valid = valid & conditions(j).validop(btd.(tp).(conditions(j).name));
                            else
                                vind = btd.findField(conditions(j).name);
                                if (vind > 0)
                                    
                                    if (strcmpi(tp, 'all'))
                                        c2 = btd.(tp).(btd.glt(vind).xField);
                                    else
                                        c2 = btd.(tp).([position '_' btd.glt(vind).xField]);
                                    end
                                    valid = valid & conditions(j).validop(btd.glt(vind).highResDerivationMethod(c2, btd.glt(vind).xData, btd.glt(vind).yData));
                                else
                                    warning ([conditions(j).name ' does not name a field in ' tp ' or a glt field']);
                                end
                            end
                        end
                        centers = centers(valid);
                    end
             end
            
            sp = interp1(xt, 1:length(xt), centers);
            
        end
            
        
        % function [bta, dm, stimnum] = behaviorTriggeredWeightedAverage (btd, tp, position, gltname, displacementAxis, varargin)
        % as with behavior triggered average, but weighted by
        % tp.(weightField)
        %optional: weightedOp -- weight is weightedOp(tp.weightField)
        function [bta, dm, stimnum] = behaviorTriggeredWeightedAverage (btd, tp, position, weightField, gltname, displacementAxis, varargin)
            % function [bta, dm] = behaviorTriggeredWeightedAverage (btd, tp, position, weightField, gltname, displacementAxis, varargin)
            weightedOp = @(x) x;
            conditions = [];
            varargin = assignApplicable(varargin);
            [dm,stimnum] = btd.behaviorTriggeredDataMatrix(tp, position, gltname, displacementAxis, 'conditions', conditions, varargin{:});
            btp = [btd.(tp)];
            wf = [btp.(weightField)];
            valid = true(size(wf));            
            for j = 1:length(conditions)
                %valid = valid & conditions(j).validop([tp.(conditions(j).name)]);
                if (isfield(btp, conditions(j).name))
                    valid = valid & conditions(j).validop([btp.(conditions(j).name)]);
                else
                    val_field = btd.behaviorTriggeredDataMatrix(tp, position, conditions(j).name, 0, varargin{:});
                    valid = valid & conditions(j).validop(val_field');
                end
            end
            wf = weightedOp(wf(valid));
            if (all(wf == mean(wf))) %#ok<CPROP>
                ww = ones(size(dm,2));
            else
%                wf = (wf-mean(wf))/std(wf); stdev removed 10/26 by MHG to
%                facilitate comparisons between data sets
                 wf = (wf-mean(wf));
                ww = repmat(wf',[1, size(dm,2)]);
            end
            dm2 = dm.*ww; 
%            ww(~isfinite(dm2)) = 0;
            dm2(~isfinite(dm2)) = 0;
            
            bta = sum(dm2)./(sum(isfinite(dm)));
            
%             dm2 = dm; dm2(~isfinite(dm2)) = 0;
%             bta = bta - sum(dm2)./sum(isfinite(dm));
%             
            stimnum = mean(stimnum);
            
        end
        
        %Takes the average of the behavior triggered data, returning the
        %average and the stretched 
        function [bta, dm] = behaviorTriggeredStretchedAverage (btd, tp, positionStart, positionEnd, gltname, npoints, varargin)
            % function [bta, dm] = behaviorTriggeredAverage (btd, tp, position, gltname, displacementAxis, varargin)
            dm = btd.behaviorTriggeredDataMatrixStretched(tp, positionStart, positionEnd, gltname, npoints, varargin{:});
            dm2 = dm; dm2(~isfinite(dm2)) = 0;
            bta = sum(dm2)./sum(isfinite(dm));
        end
        
        function [r,reb, fieldAxis] = rateVsField(btd, tp, position, gltname, fieldAxis, varargin)
            %function [r,reb, fieldAxis] = rateVsField(btd, tp, position, gltname, fieldAxis, varargin)
            %
            % r - rate of tp.position occuring vs. gltname (in Hz)
            % reb - error bar due to counting statistics
            %
            %varargin:
            % conditions: applied to tp and all
            % tpconditions: only to tp
            % allconditions: onl to all: default isrun -- if you pass
            % allconditions, make sure to include isrun if needed
            %  allconditions.name =  'isrun';allconditions.validop = @(x) logical(x);
            conditions = [];
            tpconditions = conditions;
             allconditions.name =  'isrun';allconditions.validop = @(x) logical(x);
            %conditions.name = 'hsnum';
            %conditions.validop = @(x) x == 1;
            %this will get first headsweep
            varargin = assignApplicable(varargin);
            
            yy = btd.behaviorTriggeredDataMatrix(tp, position, gltname, 0, 'conditions', [conditions tpconditions], varargin{:});
            aa = btd.behaviorTriggeredDataMatrix('all', [], gltname, 0, 'conditions', [conditions allconditions], varargin{:});
            
            dt = median(diff(btd(1).all.eti));
            
            existsAndDefault('fieldAxis', linspace(percentile(aa,0.01), percentile(aa,.99), 20));
            
            xx = binEdgesFromCenters(fieldAxis);
            
            h1 = histc(yy,xx);
            h2 = histc(aa,xx);
            rr = 1/dt*h1(1:(end-1))./h2(1:(end-1));
            reb =  1/dt*sqrt(h1(1:(end-1)))./h2(1:(end-1));
            rr= rr';
            reb = reb';
            if (nargout == 0)
                shadedErrorPlot(fieldAxis, rr, reb); ylabel([tp ' ' position ' rate (Hz)']); xlabel (gltname);
            else
                r = rr;
            end
        end
        
        function fitParams = fitRateFun(btd, tp, position, gltname, fitfun, initGuess, varargin)
            % function fitParams = fitRateFun(btd, tp, position, gltname, fitfun, initGuess, varargin)
            %
            % r - rate of tp.position occuring vs. gltname (in Hz)
            % reb - error bar due to counting statistics
            %
            %varargin:
            % conditions: applied to tp and all
            % tpconditions: only to tp
            % allconditions: only to all 
            %
            % gltname is field name or string of field names
            % fitfun = f (x, data)
            % if gltname is a cell, then data will also be a cell
            % note this behavior was revised on 9/20/2014 -- old scripts
            % may be unhappy and need to be updated
            
            conditions = [];
            tpconditions = conditions;
            allconditions.name =  'isrun';allconditions.validop = @(x) logical(x);
            %conditions.name = 'hsnum';
            %conditions.validop = @(x) x == 1;
            %this will get first headsweep
            varargin = assignApplicable(varargin);
            op = optimset('fminunc');
            op.Display = 'off';
            op.LargeScale = 'off';
            if (iscell(gltname)) 
                for j = 1:length(gltname)
                    yy{j} = btd.behaviorTriggeredDataMatrix(tp, position, gltname{j}, 0, 'conditions', [conditions tpconditions], varargin{:}); %#ok<AGROW>
                    aa{j} = btd.behaviorTriggeredDataMatrix('all', [], gltname{j}, 0, 'conditions', [conditions allconditions], varargin{:}); %#ok<AGROW>
                end
            else
                yy = btd.behaviorTriggeredDataMatrix(tp, position, gltname, 0, 'conditions', [conditions tpconditions], varargin{:});
                aa = btd.behaviorTriggeredDataMatrix('all', [], gltname, 0, 'conditions', [conditions allconditions], varargin{:});
            end
%             if (iscell(gltname))
%                 if (~iscell(fitfun) || ~iscell(initGuess))
%                     error ('all gltname, fitfun, initGuess must be cells, or none');
%                 end
%                 for j = 1:length(gltname)
%                     nargs(j) = length(initGuess{j});
%                 end
%                 
%             else
               
                myfun = @(x) -sum(log(max(1E-100,1-exp(-fitfun(x,yy))))) + sum(max(0,fitfun(x,aa)));
                
                fitParams = fminunc(myfun, initGuess, op);
%                 return;
%             end
%             function val = globalfitfun(x)
%                 argc = 0;
%                 rturn = zeros(size(yy));
%                 rall = zeros(size(aa));
%                 
%                 for k = 1:length(gltname)
%                     args = x(argc + (1:nargs(k)));
%                     argc = argc + nargs(k);
%                     rturn(:,k) = fitfun{k}(args, yy(:,k));
%                     rall(:,k) = fitfun{k}(args, aa(:,k));
%                 end
%                 val = -sum(log(max(1E-100, 1 - exp(prod(rturn, 2))))) + sum(max(0, prod(rall,2)));
%             end
%             x0 = [initGuess{:}];
%             globalfitfun(x0)
%             x1 = fminunc(@globalfitfun, x0, op);
%             globalfitfun(x1)
%             argc = 0;
%             for j = 1:length(initGuess)
%                 fitParams{j} =  x1(argc + (1:nargs(j)));
%                 argc = argc + nargs(j);
%             end
            
        end
        
        
        function btd = addConvolvedFields (btd, gltname, newname, convkernel, convkernel_dt, varargin)
            %btd = addConvolvedFields (btd, gltname, newname, convkernel, convkernel_dt, varargin)
            %NOTE: does not reverse kernel for you (e.g. convolution is in
            %normal sense)
            %NOTE: default is result of convolution is entirely causal
            % to make convolution symmetric, pass 'symmetric', true
            if ~(iscell(gltname))
                gltname = {gltname};
            end
            if ~(iscell(newname))
                newname = {newname};
            end
            existsAndDefault('convkernel_dt', []);
            symmetric = false;
            normalizeStd = false;
            existsAndDefault('scaleToSqr', false);
            varargin = assignApplicable(varargin);
            
            for j = 1:length(btd)
                if (isempty(btd(j).glt))
                    continue;
                end
                
                for k = 1:length(gltname)
                    gltind = find(strcmpi(gltname{k}, {btd(j).glt.fieldname}),1,'first');
                    if (~isempty(gltind))
                        newglt = btd(j).glt(gltind);
                        newglt.fieldname = newname{k};
                        if (isempty(convkernel_dt))
                            ck = convkernel;
                        else
                            dt = convkernel_dt(min(k, length(convkernel_dt)));
                            ck = interp1(convkernel, 1:(median(diff(newglt.xData))/dt):length(convkernel), 'linear');
                            if (scaleToSqr)
                                ck = ck*sqrt(length(convkernel)/length(ck));
                            else
                                ck = ck*length(convkernel)/length(ck);
                            end
%                            ck = interp1((1:length(convkernel))*dt, convkernel, 1::
                        end
                        if symmetric
                            newglt.yData = conv(newglt.yData, ck, 'same');
                        else
                            yy = conv(newglt.yData, ck, 'full');
                            newglt.yData = yy(1:length(newglt.yData));
                        end
%                         if (normalizeStd)
%                             newglt.yData = newglt.yData / std(newglt.yData);
%                         end
                        gltind = find(strcmpi(newglt.fieldname, {btd(j).glt.fieldname}),1,'first');
                        
                        if (isempty(gltind))
                            btd(j).glt = [btd(j).glt newglt];
                        else
                            btd(j).glt(gltind) = newglt;                           
                        end
                    end
                end
            end
            if (normalizeStd)
                for k = 1:length(newname) %#ok<UNRCH>
                    snorm = std(btd.behaviorTriggeredDataMatrix('all', '', newname{k}, 0));
                    if (isempty(snorm) || ~isfinite(snorm))
                        warning ('btd:addconvfield', 'bad data matrix, not normalizing');
                        continue;
                    end
                    for j = 1:length(btd)
                        gltind = btd(j).findField(newname{k});
                        btd(j).glt(gltind).yData = btd(j).glt(gltind).yData/snorm;
                    end
                end
            end
        end
        
        function btd = addLowpassFields (btd, gltname, sigmaTime)
            if ~(iscell(gltname))
                gltname = {gltname};
            end
            
            for j = 1:length(btd)
                if (isempty(btd(j).glt))
                    continue;
                end
                for k = 1:length(gltname)
                    gltind = find(strcmpi(gltname{k}, {btd(j).glt.fieldname}),1,'first');
                    if (~isempty(gltind))
                        newglt = btd(j).glt(gltind);
                        newglt.fieldname = [gltname{k} '_LOWPASS'];
                        sigma = sigmaTime./median(diff(newglt.xData));
                        newglt.yData = lowpass1D(newglt.yData, sigma);
                        gltind = find(strcmpi(newglt.fieldname, {btd(j).glt.fieldname}),1,'first');
                        if (isempty(gltind))
                            btd(j).glt = [btd(j).glt newglt];
                        else
                            btd(j).glt(gltind) = newglt;
                        end
                    end
                end
            end
        end
        
        
        function btd = addOperationFields (btd, gltname, newname, operation)
            % function btd = addOperationFields (btd, gltname, newname, operation)
            % if gltname and newname are both single names of fields
            % or if gltname and newname are cells of the same length
            % we create fields s.t. newname{j}.yData = operation(gltname{j}.yData)
            % if gltname is a cell and newname is a single field, then
            % newname = operation({gltname{1}.yData, gltname{2}.yData,
            % ...})
            % in this case, operation should take a cell of data
            % TODO: REVISE TO BE CORRECT FUNCTION FOR several glt, 1
            % newname
            if ~(iscell(gltname))
                gltname = {gltname};
            end
            if ~(iscell(newname))
                newname = {newname};
            end
            
            for j = 1:length(btd)
                if (isempty(btd(j).glt))
                    continue;
                end
                
                if (length(newname) == length(gltname))
                    for k = 1:length(gltname)
                        gltind = find(strcmpi(gltname{k}, {btd(j).glt.fieldname}),1,'first');
                        if (~isempty(gltind))
                            newglt = btd(j).glt(gltind);
                            newglt.fieldname = newname{k};
                            newglt.yData = operation(newglt.yData);                            
                        else
                            warning ('btd:addoperation', [gltname{k} ' not found']);
                            continue;
                        end
                        gltind = find(strcmpi(newglt.fieldname, {btd(j).glt.fieldname}),1,'first');
                        if (isempty(gltind))
                            btd(j).glt = [btd(j).glt newglt];
                        else
                            btd(j).glt(gltind) = newglt;
                        end
                    end
                else
                    gltind = btd(j).findField(gltname{1});
                    if (gltind <= 0)
                        error ([gltname{1} ' not found']);
                    end
                    newglt = btd(j).glt(gltind);
                    
                    newglt.fieldname = newname{1};
                    yData{1} = newglt.yData;
                    
                    for k = 2:length(gltname)
                        gltind = btd(j).findField(gltname{k});
                        if (gltind <= 0)
                            error ([gltname{k} ' not found']);
                        end
                        yData{k} = interp1(btd(j).glt(gltind).xData, btd(j).glt(gltind).yData, newglt.xData, 'linear', 'extrap'); %#ok<AGROW>
                    end
                    newglt.yData = operation(yData);
                    gltind = find(strcmpi(newglt.fieldname, {btd(j).glt.fieldname}),1,'first');
                    if (isempty(gltind))
                        btd(j).glt = [btd(j).glt newglt];
                    else
                        btd(j).glt(gltind) = newglt;
                    end
                end
                
                
            end
            
        end
        
        function  [datamatrix, stimnum] = baseConditionedDataMatrix(btd, gltname, conditions, displacementAxis,varargin)
        %function  [datamatrix, stimnum] = baseConditionedDataMatrix(btd, gltname, conditions, displacementAxis,varargin)
        %
        % datamatrix of gltname centered on points chosen from gltname.xData satisfying
        % conditions(j).validop(conditions(j).name at point) is true for all j
        % btd < BehaviorTriggeredData
        % gltname < name of a glt in btd
        % conditions < array of conditions: 
        % conditions(j).name is name of a glt in btd
        % conditions(j).validop is operation that takes yData from glt and
        % returns true/false
        % displacementAxis - axis to take average relative to chosen center
        % points
        % optional arg: maxDmElems, [1E7]
        % to avoid taking forever, center points are decimated to keep
        % datamatrix under maxDmElems per btd per gltname

            if (length(btd) > 1)
                [datamatrix, stimnum] =  baseConditionedDataMatrix (btd(1), gltname, conditions, displacementAxis, varargin{:});
                for j = 2:length(btd)
                    if (isempty(btd(j).glt))
                        continue;
                    end
                    [dm,sn] =  baseConditionedDataMatrix (btd(j), gltname, conditions, displacementAxis, varargin{:});
                    datamatrix = [datamatrix;dm]; %#ok<AGROW>
                    stimnum = [stimnum;sn];%#ok<AGROW>
                end
                return;
            end
            
            if (iscell(gltname))
                datamatrix =  baseConditionedDataMatrix (btd, gltname{1}, conditions, displacementAxis, varargin{:});
                stimnum = ones(size(datamatrix));
                for k = 2:length(gltname)
                    datamatrix = [datamatrix  baseConditionedDataMatrix(btd, gltname{k}, conditions, displacementAxis, varargin{:})]; %#ok<AGROW>
                    stimnum = [stimnum k*ones(size(datamatrix))]; %#ok<AGROW>
                end
                return;
            end
            maxDmElems = 1E7;
            varargin = assignApplicable(varargin);
            maxCenterPoints = ceil(maxDmElems/length(displacementAxis));
            
            %Make a datamatrix from the glt
            ind = btd.findField(gltname);
            if (ind <= 0)
                datamatrix = [];
                stimnum = [];
                return;
            end
            xd = btd.glt(ind).xData;
            valid = true(size(xd));
            
            try
                for j = 1:length(conditions)
                    cind = btd.findField(conditions(j).name);
                    yd = btd.glt(cind).highResDerivationMethod(xd, btd.glt(cind).xData, btd.glt(cind).yData);
                    valid = valid & conditions(cind).validop(yd);
                end
            catch me
                disp(me.getReport());
                datamatrix = [];
                stimnum = [];
                return;
            end
            
            centers = xd(valid);
            if (length(centers) > maxCenterPoints)
                centers = centers(randperm(length(centers), maxCenterPoints));
            end
            datamatrix = btd.glt(ind).makeDataMatrix(centers, displacementAxis, varargin{:});
            stimnum = ones(size(datamatrix));
        end
    
        function  [bca, dm, stimnum] = baseConditionedAverage(btd, gltname, conditions, displacementAxis,varargin)
        %function  [bca, datamatrix, stimnum] = baseConditionedAverage(btd, gltname, conditions, displacementAxis,varargin)
        %
        % average of gltname centered on points chosen from gltname.xData satisfying
        % conditions(j).validop(conditions(j).name at point) is true for all j
            [dm, stimnum] = baseConditionedDataMatrix(btd, gltname, conditions, displacementAxis,varargin{:});
            dm2 = dm; dm2(~isfinite(dm2)) = 0;
            bca = sum(dm2)./sum(isfinite(dm));
            stimnum = mean(stimnum);
        end
        
        function [cm, cm_stim, stimnum, dm] = behaviorTriggeredCovarianceMatrix (btd, tp, position, gltname, displacementAxis, varargin)
            
            if ~(iscell(gltname))
                gltname = {gltname};
            end
            
            dm = behaviorTriggeredDataMatrix(btd, tp, position, gltname, displacementAxis, varargin{:});
            
            dm2 = dm; dm2(~isfinite(dm)) = 0;
            cm = dm2'*dm2 / size(dm2,1);
            cm_stim = zeros(size(cm));
            
            stimnum = zeros(size(dm));
            for j = 1:length(gltname)
                dl = length(displacementAxis);
                stimnum(:, (j-1)*dl + (1:dl)) = j;
            end
            
            %estimate covariance matrices of signal
            btdw = zeros(size(btd));
            for j = 1:length(btd)
                if (isempty(btd(j).glt))
                    continue;
                end
                btdw(j) = length(btd(j).all.eti);
                for k = 1:length(gltname)
                    for m = 1:length(gltname)
                        glt1 = btd(j).glt(find(strcmpi([gltname{k}], {btd(j).glt.fieldname}),1,'first'));
                        glt2 = btd(j).glt(find(strcmpi([gltname{m}], {btd(j).glt.fieldname}),1,'first'));
                        
                        yd2 = interp1(glt2.xData, glt2.yData, glt1.xData, 'linear', 0);
                        dt = median(diff(glt1.xData));
                        maxlag = ceil(2*max(abs(displacementAxis))/dt);
                        [xc, lags] = xcorr(glt1.yData, yd2, maxlag, 'unbiased');
                        off1 = (k-1)*length(displacementAxis);
                        off2 = (m-1)*length(displacementAxis);
                        for n = 1:length(displacementAxis)
                            cm_stim (off1 + n, off2 + (1:length(displacementAxis))) = cm_stim (off1 + n, off2 + (1:length(displacementAxis))) + btdw(j)*interp1(lags*dt, xc, displacementAxis - displacementAxis(n));
                        end
                        
                    end
                end
            end
            cm_stim = cm_stim / sum(btdw);
        end
        
        function [k0,k1,k2] = wienerKernels(btd, gltname, dqname, displacementAxis, varargin)
        %function [k0,k1,k2] = wienerKernels(gltname, dqname, displacementAxis, varargin)

            da = sort(unique(displacementAxis), 'ascend');
            dt = min(diff(da));
            dtf = dt/2;
            aa = [btd.all];
            
            for j = 1:length(aa)
                y{j} = aa(j).(dqname);                
                eti{j} = aa(j).eti;
                dtf = max(dtf, median(diff(aa(j).eti)));
            end
            k0 = mean([y{:}]);
            if (nargout < 2)
                return;
            end
             maxlag = ceil((max(da) - min(da))/(2*dtf)) + 1;
            c = zeros(length(aa), 2*maxlag + 1);
            ns = c;
            v = zeros(1,length(aa));
            for j = 1:length(aa)
                ind = btd(j).findField(gltname);
                if (ind <= 0)
                    warning ('btd:wienerKernels', [gltname ' is not a glt']);
                    k1 = []; k2 = [];
                    return;
                end
                glt = btd(j).glt(ind); %#ok<PROP>

                midtime = median(da);
                txf = max(min(eti{j} - midtime/2), min(glt.xData + midtime/2)):dtf:min(max(eti{j} -midtime/2),max(glt.xData + midtime/2)); %#ok<PROP>
                [xx,~,sumy, numx] = meanyvsxFast(eti{j} - midtime/2, y{j}-k0, txf);
                
                u = GlobalLookupTable.averageInPrecedingBin(xx, glt.xData + midtime/2, glt.yData);	
               
            
                [c(j,:),lags] = xcorr(sumy,u,maxlag,'none');
                ns(j,:) = xcorr(numx, ones(size(u)), maxlag, 'none');
                v(j) = var(glt.yData);
            end
            cc = sum(c, 1)./sum(ns, 1);
            n2 = sum(ns, 2);
            vu = v*n2/sum(n2);

            k1 = interp1(lags*dtf + midtime, cc/vu, displacementAxis, 'linear');
            
            %only do the convolution in the causal direction -- we should
            %really clean this code up
            da = displacementAxis(displacementAxis >= 0);
            kcausal = k1(displacementAxis >= 0);
            %y1 = conv(yd, cc, 'full'); y1 = y1(1:length(yd)); sppred = interp1(xd, y1, btd(1).all.eti);[~,my2] = meanyvsxFast(btd(1).all.eti,sppred, 0:0.1:1200);[c,lags] = xcorr(my, my2, 300); plot (lags/10, c)
            
            %vu2
            if (nargout < 3)
                return;
            end
            
            warning ('k2 not implemented yet'); 
            
            
            return;
            
        end

        function [mi, tcent, bit_rate, tstart] = miFieldSignalVsTime (btd, tp, field, position, gltname, taxis, tbinEdges, varargin)
        %function [mi, tcent, bit_rate, tstart] = miFieldSignalVsTime (btd, tp, field, position, gltname, conditions, taxis, tbinEdges, varargin)
            conditions = [];
            operation = @(x) double(x);
            method = 'Shannon_Edgeworth';
            varargin = assignApplicable(varargin);
            aa = [btd.(tp)];
            fd = operation([aa.(field)]);
            valid = true(size(fd));
            for j = 1:length(conditions)
                valid = valid & conditions(j).validop([aa.(conditions(j).name)]);
            end
            fd = fd(:,valid);
            dm = btd.behaviorTriggeredDataMatrix(tp, position, gltname, taxis, 'conditions', conditions, varargin{:});
  
            if (size(tbinEdges,1) == 1)
               if (size(tbinEdges,2) == 1) %it's a bin width rather than edges
                    bw = tbinEdges;
                    tbinEdges = taxis(taxis < max(taxis - bw));
                    tbinEdges(2,:) = tbinEdges(1,:) + bw;
               else
                   tbinEdges(2,:) = [tbinEdges(2:end) tbinEdges(end)+eps];
               end
           end
           
          
           mi = zeros(1, size(tbinEdges,2)-1);
           tcent = mi;
           co = IShannon_HShannon_initialization(1, {'member_name', method});
           for j = 1:(size(tbinEdges,2)-1)
               inds = taxis >= tbinEdges(1,j) & taxis < tbinEdges(2,j);
               tcent(j) = mean(taxis(inds));
               ds = [size(fd,1) nnz(inds)];
               yy = [fd' dm(:,inds)];
               yy = yy(all(isfinite(yy),2), :); %#ok<CPROP> all is also a property; function meant here
               if (isempty(yy))
                   continue;
               end
               mi(j) = IShannon_HShannon_estimation(yy', ds, co);
           end
           bit_rate = (mi/log(2)) ./ diff(tbinEdges(:,1:end-1));
           tstart = tbinEdges(1,1:end-1);
               
        end
        
        function [mi, tcent, bit_rate, tstart] = miVsTime (btd, items, taxis, tbinEdges, varargin)
            %function [mi, tcent] = miVsTime (btd, items, taxis, tbinEdges, varargin)
            %finds the mutual information between signals defined by
            %items.gltname and states defined by different items
            %see static function miDiscreteVsTime
            %
            %example
            % firsths.name = 'hsnum';
            % firsths.validop = @(x) x == 1;
            % items(1).tp = 'rejhs';
            % items(1).position = 'acchs';
            % items(1).gltname = 'led1ValDiff';
            % items(1).conditions = firsths;
            % items(2) = items(1);
            % items(2).tp = 'acchs';
            %
            % binWidth = 0.5;
            % [mi,tcent] = btd.miVsTime(items, -4:0.05:4, binWidth);
            % mi_rate = mi/(log(2)*binWidth); %information rate about hs rejection/acceptance in bits/sec about tcent
            
            if (~exist('items', 'var'))
                items.tp = '';
                items.position = '';
                items.gltname = '';
                items.conditions.name = '';
                items.conditions.validop = '';
                mi = items;
                return;
            end
            if ~isfield(items, 'conditions')
                [items.conditions] = deal([]);
            end
            for j = 1:length(items)
                %(btd, tp, position, gltname, displacementAxis, varargin)
                dm{j} = btd.behaviorTriggeredDataMatrix(items(j).tp, items(j).position, items(j).gltname, taxis, 'conditions', items(j).conditions, varargin{:});
            end
            [mi, tcent, bit_rate, tstart]  = BehaviorTriggeredData.miDiscreteVsTime(taxis, dm, tbinEdges, varargin{:});
          
        end
       
        function [dkl12, tcent, bit_rate, tstart] = klDivergenceVsTime (btd, item1, item2, taxis, tbinEdges, varargin)
            %function [mi, tcent] = miVsTime (btd, items, taxis, tbinEdges, varargin)
            %finds the kl divergence dkl(item1||item2) information between signals defined by
            %item1 and item2 
            %see static function klDivergenceBetweenDM
            %
           
            if (~exist('item1', 'var'))
                items.tp = '';
                items.position = '';
                items.gltname = '';
                items.conditions.name = '';
                items.conditions.validop = '';
                dkl12 = items;
                return;
            end
            if ~isfield(item1, 'conditions')
                item1.conditions = [];
            end
            if ~isfield(item2, 'conditions')
                item2.conditions = [];
            end
            items = [item1, item2];
            for j = 1:length(items)
                %(btd, tp, position, gltname, displacementAxis, varargin)
                dm{j} = btd.behaviorTriggeredDataMatrix(items(j).tp, items(j).position, taxis, 'conditions', items(j).conditions, varargin{:});
            end
             % function [dkl12, tcent] = klDivergenceBetweenDM (taxis, dm1, dm2, tbinEdges, varargin)
            [dkl12, tcent, bit_rate, tstart]  = BehaviorTriggeredData.klDivergenceBetweenDM(taxis, dm1, dm2, tbinEdges, varargin{:});
            
        end
    
        function [fom_true, fom_shifted, u, s, vals, ue, se] = figureOfMerit (btd, tp, position, gltname, displacementAxis, varargin)
            %[fom_true, fom_shifted, u, s, vals, ue, se] = figureOfMerit (btd, tp, position, gltname, displacementAxis, varargin)
            %arguments same as for behaviorTriggeredAverage
            %additional arguments
            %shiftPad - added to minimum shift to exclude slightly larger
            %region around 0 shift
            %deltaT - amount to shift at one time - default duration of
            %displacementAxis
          
            shiftPad = 0;
            deltaT = [];
            varargin = assignApplicable(varargin);
            
            [~,~,bn] = btd.behaviorTriggeredDataMatrix(tp, position, gltname, displacementAxis, varargin{:});
            nturns = zeros([1 length(btd)]);
            for j = 1:length(btd)
                nturns(j) = nnz(bn == j);
            end
            [ue,se] = btd.gatherGlt(gltname).theoreticalSumSqStats(max(displacementAxis)-min(displacementAxis), length(displacementAxis), nturns);

            trange = 0;
            for j = 1:length(btd)
                centers = btd(j).getBTDMCenters(tp, position, gltname, varargin{:});
                trange = max(trange, max(centers)-min(centers));
            end
            if (isempty(deltaT))
                deltaT = max(displacementAxis) - min(displacementAxis);
            end
            nshifts = floor((trange-2*shiftPad)/deltaT)-1;
            
            vals = zeros([1 nshifts]);
            for j = 1:nshifts
                vals(j) = sum(btd.behaviorTriggeredAverage(tp, position, gltname, displacementAxis, 'shift', j*deltaT+shiftPad, varargin{:}).^2);
            end
            u = mean(vals);
            s = std(vals);

            fom_true = (sum(btd.behaviorTriggeredAverage(tp, position, gltname, displacementAxis, varargin{:}).^2) - ue)/se;
            fom_shifted = (vals - ue)/se;
            
        end
    
        function [rate, ttime, rtime] = simulateTurnTrain (btd, xname, ratefun, taxis, nlarvae, nreps)
            % function [rate, ttime, rtime] = simulateTurnTrain (btd, xname, ratefun, taxis, nlarvae)
            % rate = ratefun({xdata}, t) -- {xdata} is a cell of values at
            % times specified in taxis; i.e. xdata{j} = xname{j}(taxis)
            % ttime is time of turns (given nlarvae and poisson statistics)
            % rtime is time of runs (currently just taxis repeated nlarvae
            % times)
            % simulation of turns is repeated nreps times
            if (length(btd) > 1)
                for j = 1:length(btd)
                    [ttime{j}, rtime{j}] = simulateTurnTrain (btd(j), xname, ratefun, taxis, nlarvae, nreps); %#ok<AGROW>
                end
                return;
            end
            if (~iscell(xname))
                xname = {xname};
            end
            %gather xdata
            for j = 1:length(xname)
                ind = btd.findField(xname{j});
                if (ind <= 0 && j == 1)
                    warning ([xname{j} ': field not found']);
                    rtime = [];
                    ttime = [];
                    return;
                end
                xf = btd.glt(ind);
                if (~strcmpi(xf.xField, 'eti'))
                    warning (['eti is not xField for ' xf.fieldname ' results possibly incorrect']);
                end
                existsAndDefault('taxis', xf.xData);
                xdata{j} = xf.derivationMethod(taxis, xf.xData, xf.yData); %#ok<AGROW>
            end
            
            rate = ratefun(xdata, taxis);
            if (nargout < 2)
                return;
            end
            existsAndDefault('nlarvae', 1);
            existsAndDefault('nreps', 1);
            
            dt = diff(taxis); dt = [dt(1) dt];
            meanturns = nlarvae*rate.*dt;
%             plot (taxis, meanturns); pause;
            for j = 1:nreps
                nturns = poissrnd(meanturns);
                nt = nturns(nturns > 0);
                tt = taxis(nturns > 0);
                ttime{j} = []; %#ok<AGROW>
                for k = 1:max(nt)
                    ttime{j} = [ttime{j} tt(nt >= k)];
                end
                ttime{j} = sort(ttime{j}(:)); %#ok<AGROW>
                
            end
            if (nreps == 1)
                ttime = ttime{1};
            end
            rtime = repmat(taxis, [nlarvae 1]); rtime = rtime(:);
            
            
        end
        
         problemDescription = simulateVarSwitch (btd, opstruct, oldpd, timeField, trange, staticParams, temporalParams, nlarvae);    
    end
end


