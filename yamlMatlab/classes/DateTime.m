classdef DateTime
    % This class enclapsulates date time value but behaves in a very
    % similar way as typical array of timestamps in datenum format
  %======================================================================
    %{
		Copyright (c) 2011
		This program is a result of a joined cooperation of Energocentrum
		PLUS, s.r.o. and Czech Technical University (CTU) in Prague.
        The program is maintained by Energocentrum PLUS, s.r.o. and 
        licensed under the terms of MIT license. Full text of the license 
        is included in the program release.  		
		
        Author(s): 
		Jiri Cigler, Dept. of Control Engineering, CTU Prague & Automatic Control Laboratory, ETH Zurich		
		Jan  Siroky, Energocentrum PLUS s.r.o.		
		
        Implementation and Revisions:

        Auth  Date        Description of change
        ----  ---------   -------------------------------------------------
        jc    01-Mar-11   First implementation
    %}
    %======================================================================

    properties
        serialDate
    end
    methods
        function this = DateTime(varargin)
            if numel(varargin)==1 && isa(varargin{1},'java.util.Date')
                    this.serialDate=-java.util.Date().getTimezoneOffset()/(60*24)+datenum(char(varargin{1}.toGMTString));
            else
                this.serialDate=datenum(varargin{:});
            end
        end
        function this = plus(this,val)
            o =@plus;
            this = doFun(this,o,val);
        end
        function this = minus(this,val)
            o =@minus;
            this = doFun(this,o,val);
        end
        function this = times(this,val)
            o =@times;
            this = doFun(this,o,val);
        end
        
        function this = mtimes(this,val)
            o =@mtimes;
            this = doFun(this,o,val);
        end
        
        function this = mrdivide(this,val)
            o =@mrdivide;
            this = doFun(this,o,val);
        end
        
        function this = rdivide(this,val)
            o =@rdivide;
            this = doFun(this,o,val);
        end
        
        
        
        function this = horzcat(this,varargin)
            %this.serialDate = [this.serialDate, n.serialDate];
            for i=1:numel(varargin)
                this.serialDate = [this.serialDate, varargin{i}.serialDate];
            end
        end
        
        function this = vertcat(this,varargin)
            for i=1:numel(varargin)
                this.serialDate = [this.serialDate; varargin{i}.serialDate];
            end
        end
        
        
        function this = ctranspose(this)
            this.serialDate = this.serialDate';
        end
        
        function this = transpose(this)
            this.serialDate = this.serialDate';
        end
        function  disp(this)
            disp([this.serialDate])
        end
        function out = double(this)
            out = this.serialDate;
        end
        function out = length(this)
            out = length(this.serialDate);
        end
        
        function out = size(this,varargin)
            out = size(this.serialDate,varargin{:});
        end
        
        function out = numel(this)
            out = numel(this.serialDate);
        end
        function out = isreal(this)
            out = isreal(this.serialDate);
        end
        function out = isnan(this)
            out = isnan(this.serialDate);
        end
        function out = isfinite(this)
            out = isfinite(this.serialDate);
        end
        
        function out = le(this,B)
            if isa(B,'DateTime')
                out = le(this.serialDate,B.serialDate);
            else
                out = le(this.serialDate,B);
            end
        end
        
        function out = lt(this,B)
            fun=@lt;
            if isa(B,'DateTime')
                out = fun(this.serialDate,B.serialDate);
            else
                out = fun(this.serialDate,B);
            end
        end
        function out = gt(this,B)
            fun=@gt;
            if isa(B,'DateTime')
                out = fun(this.serialDate,B.serialDate);
            else
                out = fun(this.serialDate,B);
            end
        end
        function out = eq(this,B)
            fun=@eq;
            if isa(B,'DateTime')
                out = fun(this.serialDate,B.serialDate);
            else
                out = fun(this.serialDate,B);
            end
        end
        function out = diff(this)
            out = diff(this.serialDate);
        end
        
        function out = norm(this,varargin)
            out = norm(this.serialDate,varargin{:});
        end
        
        function [this k] = sort(this,varargin)
            [this.serialDate k] = sort(this.serialDate,varargin{:});
        end
        
        function this = subsref(this,S)
            if isa(S.subs{1},'DateTime')
                S.subs{1}=double(S.subs{1});
            end
                
            this.serialDate =  subsref(this.serialDate,S);
        end
        
        function idx = subsindex(this)
            idx = double(this)-1;
        end
        
        function endidx = end(this,k,n)  
            if size(this.serialDate,1)==1 || size(this.serialDate,2)==1
                endidx=numel(this.serialDate);
            else
                endidx = size(this.serialDate,k);
            end
        end
        
        function this = subsasgn(this, S, B)
            if not(isa(B,'DateTime'))
                B=DateTime(B);
            end
            
            this.serialDate =subsasgn(this.serialDate, S, B);
        end
        
        function res = bsxfun(fun,A,B)
            res = fun(A,B);
        end
        
        function out =superiorfloat (x,y,xi)
            if isa(x,'DateTime') && isa(xi,'DateTime')
                out = superiorfloat(x.serialDate,y,xi.serialDate);
            elseif isa(x,'DateTime') && not(isa(xi,'DateTime'))
                out = superiorfloat(x.serialDate,y,xi);
            elseif not(isa(x,'DateTime')) && isa(xi,'DateTime')
                out = superiorfloat(x,y,xi.serialDate);
            else
                out = superiorfloat(x,y,xi);
            end
        end
        
        function this = floor(this)
            this.serialDate = floor(this.serialDate);
        end
        function this = max(this,varargin)
            this.serialDate = max(this.serialDate,varargin{:});
        end
        function this = min(this,varargin)
            this.serialDate = min(this.serialDate,varargin{:});
        end
        function out = datestr(this)
            out = datestr(this.serialDate);
        end
    end
    
    methods (Access = private)
        function this = doFun (this,o, val)
            if isa(val,'DateTime') && isa(this,'DateTime')
                this.serialDate=o(this.serialDate, val.serialDate);
            elseif isa(val,'DateTime') && not(isa(this,'DateTime'))
                val.serialDate=o(this, val.serialDate);
                this = val;
            elseif not(isa(val,'DateTime')) && (isa(this,'DateTime'))
                this.serialDate=o(this.serialDate, val);
            else
                this.serialDate=DateTime(o(this, val));
            end
        end
        
        
    end
    
end
