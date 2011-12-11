%VSIZE	disect a variable and show the size of its members
%	VSIZE shows the size of the contents of a ML variable
%	and all of its containers/members
%	see also: who, whos, ndims, size, numel, length
%
%SYNTAX
%-------------------------------------------------------------------------------
%		    VSIZE(VAR1,VAR2,...,VARn)
%		    VSIZE(OPT,VAR1,VAR2,...,VARn)
%		P = VSIZE(...)
%
%INPUT
%-------------------------------------------------------------------------------
% VARx	:	a valid, named or an unnamed ML variable
% OPT		processing mode (MUST be first arg!)
% --------------------------------------------------
% -c	:	show container classes     only
% -r	:	show top-level variable(s) only
%
%OUTPUT
%-------------------------------------------------------------------------------
% P	:	a cell matrix with all information gathered during processing
% col 1	:	field descriptor
% col 2	:	size
% col 3	:	class
% col 4	:	member descriptor			NDIM:SIZE:CLASS
% col 5	:	depth  level of container(s)/member(s)
% col 6	:	number of input argument(s)
% col 7	:	number of container(s)			cell|struct
% col 8	:	number of member(s)			all other objects
%
%EXAMPLE
%-------------------------------------------------------------------------------
%		clear v;
%		v(1).a{1}=sparse(magic(3)+2i*magic(3));
%		v(2).a{2}={struct('FA',{'a','bb'},'FB',{magic(5),{}})};
%		v(2).b{2}=@(x) sind(x);
%		e=vsize(v,v(1).a{1});
% ----------------------
%       1346       1346*   v = 2:1x2:struct
% CELL -----        256    v[].a = 2:1x1:cell
%       1150        196-   v[].a{} = 2:3x3:double.sparse.complex
% CELL -----        698    v[].a = 2:1x2:cell
%       1150          0-   v[].a{} = 2:0x0:double
% CELL -----        634    v[].a{} = 2:1x1:cell
% STRUCT ---        574    v[].a{}{} = 2:1x2:struct
%       1148          2-   v[].a{}{}[].FA = 2:1x1:char
%       1144          4-   v[].a{}{}[].FA = 2:1x2:char
%        944        200-   v[].a{}{}[].FB = 2:5x5:double
% CELL -----          0    v[].a{}{}[].FB = 2:0x0:cell
%        944          0-   v[].b = 2:0x0:double
% CELL -----         80    v[].b = 2:1x2:cell
%        944          0-   v[].b{} = 2:0x0:double
%        928         16-   v[].b{} = 2:1x1:function_handle
% ----------------------
%        196        196*   * = 2:3x3:double.sparse.complex
%
% col 1:	top-level size - member size [bytes]
%		- container sizes are NOT subtracted
%		- the last number shows the overhead of the containers
%			argument 1:	928
%			argument 2:	has no containers|members
% col 2:	size of the member [bytes]
%		###* = top-level
%		###- = member
% col 3:	field descriptor = member dimension:size:class
%		- * = unnamed input
%
%		vsize('-c',v,1i);
% ----------------------
%       1346       1346*   v = 2:1x2:struct
% CELL -----        256    v[].a = 2:1x1:cell
% CELL -----        698    v[].a = 2:1x2:cell
% CELL -----        634    v[].a{} = 2:1x1:cell
% STRUCT ---        574    v[].a{}{} = 2:1x2:struct
% CELL -----          0    v[].a{}{}[].FB = 2:0x0:cell
% CELL -----         80    v[].b = 2:1x2:cell
% ----------------------
%         16         16*   * = 2:1x1:double.complex

% created:
%	us	08-Aug-2006
% modified:
%	us	25-Feb-2009 00:40:39
% localid:	us@USZ|ws-nos-36362|x86|Windows XP|7.7.0.471.R2008b

function	[p,par]=vsize(varargin)

		earg={
			'-r'	% top-level only
			'-c'	% containers only
		};
		cc={
			'cell'
			'struct'
		};

	if	nargout
		p=[];
	end
	if	nargin < 1
		help(mfilename);
		return;
	end

% initialize engine
		[p,par]=set_par(-1,{earg,cc},nargin,[],varargin{:});
	for	i=par.nb:nargin
% tokenize input arguments
		inam=inputname(i);
		[p,par]=run_vsize(0,inam,p,par,varargin{i});
	end
		p=par.r;
	if	~nargout
		clear p;
	end
end
%-------------------------------------------------------------------------------
function	[p,par]=run_vsize(isrec,inam,p,par,varargin)

		[p,par]=set_par(isrec,inam,p,par,varargin{:});
	if	isempty(p)
		return;
	end
		[p,par]=get_tok(p,par);
end
%-------------------------------------------------------------------------------
function	[p,par]=set_par(isrec,inam,p,par,varargin)
% create common parameter list

			par.darg=4;
			par.narg=0;
			par.hasvar=false;
			par.hasanomaly=false;
			par.fnam=[];

	switch	isrec
% recursion
	case	1
			par.cflg=true;
			par.rflg=true;
			par.level=varargin{1};
			par.vnam=varargin{2};
			par.arg=varargin(3:end);
			par.narg=numel(par.arg);
			par.inam=par.vnam;
		if	par.level >= par.rlim
			disp(sprintf('VSIZE> recursion limit reached %5d',par.rlim));
			p=[];
		end
			return;

% input argument
	case	0
			par.cflg=false;
			par.sflg=0;
			par.noname='*';
			par.level=0;
			par.b=0;
			par.inam=par.noname;
		if	~isempty(inam)
			par.inam=inam;
		end
			par.vnam=par.inam;
			par.arg=varargin;
			par.narg=numel(par.arg);
			par.na=par.na+1;
			par.ne=1;
			par.nm=0;

% initialize engine
	case	-1
			par.opt=inam{1};
			par.cc=inam{2};

			narg=p;
			p=[];
			p.narg=narg;
			par.cflg=false;
			par.rflg=true;
			par.Cflg=false;
			par.nb=1;
		if	ischar(varargin{1})		&&...
			varargin{1}(1)=='-'
		if	strcmp(par.opt{1},varargin{1})	% -r
			par.rflg=false;
			par.nb=2;
		end
		if	strcmp(par.opt{2},varargin{1})	% -c
			par.Cflg=true;
			par.nb=2;
		end
		end
			par.sflg=0;
			par.rlim=get(0,'recursionlimit');
			par.argn=nargin-par.darg;
			par.na=0;
			par.ne=0;
			par.nm=0;
			par.r={};

% - output format
			nf=3-par.rflg;
			fmt={
				'%3d>%5d/%5d/%5d: %10s %10d%s   %s = %s'	1	45
				'%% %10s %10d%s   %s = %s'				5	22
				'%10d%s   %s = %s'				6	12
			};
			par.fmt=fmt{nf,1};
			par.fmtb=fmt{nf,2};
			par.d=['% ',repmat('-',1,fmt{nf,3})];

% - get current <whos> entries
%	7.2.0.232	(R2006a)
%	7.3.0.32269	(R2006b)
%	7.7.0.471	(R2008b)
			wn={
				'name'
				'size'
				'bytes'
				'class'
%				'global'
%				'sparse'
%				'complex'
				'nesting'
%				'persistent'
			};
			w=whos('nan');
			wf=fieldnames(w);
			par.wf=wf(~ismember(wf,wn));
			par.anomaly='2:0x0:anomaly (struct)';
	end
end
%-------------------------------------------------------------------------------
function	[p,par]=get_tok(p,par)
% tokenize input arguments

	for	n=1:par.narg
			arg=par.arg{n};
			[p,par]=get_arg(p,par,arg);
			[p,par]=get_fld(p,par,arg);
			[p,par]=get_ent(p,par,arg);
	end
end
%-------------------------------------------------------------------------------
function	[p,par]=get_arg(p,par,arg)

			par.hasvar=false;
	if	isstruct(arg)
			par.fnam=fieldnames(arg);
			par.sflg=par.sflg+1;
		if	par.cflg		&&...
			par.sflg ~= 1
			par.hasvar=true;
			par.sflg=0;
		end
		if	par.sflg==1
			par.hasvar=false;
			par.sflg=par.sflg+1;
		end
	else
			par.sflg=0;
	end

			a=whos('arg');
		if	~par.level
			par.b=a.bytes;
		end

	if	par.hasvar
			p.nv=numel(arg);
			p.nf=numel(par.fnam);
			p.n=p.nv*p.nf;
			p.as=size(arg);
			p.fr=cell(p.n,1);
			p.fn=cell(p.n,4);
			ix=0;
		for	j=1:p.nf
		for	i=1:p.nv
			ix=ix+1;
			p.fr(ix,1)=par.fnam(j);
			p.fn(ix,1)=par.fnam(j);
		end
		end
	else
			p.nv=1;
			p.nf=numel(par.fnam);
			p.n=1;
			p.as=size(arg);
			p.fn=cell(1,4);
			p.fr={[]};
	end
end
%-------------------------------------------------------------------------------
function	[p,par]=get_fld(p,par,arg)

	for	i=1:p.n
			vx=rem(i-1,p.nv)+1;
	if	par.hasvar
			cf=p.fn{i,1};
		try
			t=arg(vx).(cf);
		catch	%#ok
			disp('VSIZE> unexpected assignment error!');
			keyboard
		end
			p.fn{i}=sprintf('%s.%s',par.inam,p.fn{i});
	else
			t=arg;
			p.fn{i}=sprintf('%s',par.inam);
	end

			w=whos('t');
			c=w.class;
		for	j=1:numel(par.wf)
		if	w.(par.wf{j})
			c=[c,'.',par.wf{j}];		%#ok
		end
		end
			spc=sprintf('%-1dx',w.size);
			spc=sprintf('%-1d:%s',ndims(t),spc);
			spc=sprintf('%s:%s',spc(1:end-1),c);
			p.fn{i,2}=w.bytes;
			p.fn{i,3}=w.class;
			p.fn{i,4}=spc;

		if	isnumeric(p.as)			&&...
			~sum(p.as)
			par.hasanomaly=true;
			p.fn{i,3}=class(arg);
		end
	end
		if	isempty(p)			||...
			~p.n
			par.hasanomaly=true;
			p.n=numel(par.fnam);
			p.fr={[]};
		for	i=1:p.n
			p.fn(i,:)={[par.inam,'.',par.fnam{i}],0,'anomaly',par.anomaly};
		end
		end
end
%-------------------------------------------------------------------------------
function	[p,par]=get_ent(p,par,arg)

	for	i=1:p.n
			sflg=true;
		if	par.level
% a member
		if	~any(ismember(par.cc,p.fn{i,3}))
			par.b=par.b-p.fn{i,2};
			par.db=sprintf('%-1d',par.b);
			par.dlevel='-';
			par.nm=par.nm+1;
		if	par.Cflg
			sflg=false;
		end
		else
% a container
			cc=p.fn{i,3};
			par.db='container';
			par.db=sprintf('%-9.9s',upper(cc));
			par.db(numel(cc)+2:10)='-';
			par.dlevel=' ';
			par.ne=par.ne+1;
			par.nm=0;
		end
		else
% top-level variable
			par.db=sprintf('%-1d',par.b);
			par.dlevel='*';
		if	par.rflg
			disp(par.d);
		end
		end

% update/display current structure
			par=show_ent(sflg,p,par,i,par.level,par.na,par.ne,par.nm,par.db,p.fn{i,2},par.dlevel,p.fn{i,1},p.fn{i,4});
			par.hasanomaly=false;
			[p,par]=get_rec(p,par,arg,i);
	end
end
%-------------------------------------------------------------------------------
function	[p,par]=get_rec(p,par,arg,ix)
% resolve recursion

	if	par.rflg				&&...
		p.fn{ix,2}

			par.level=par.level+1;
			vx=rem(ix-1,p.nv)+1;
	switch	p.fn{ix,3}
% - struct
	case	'struct'
			par.sflg=-1;
			sb='';
		if	~isempty(p.fr{ix})
			val=arg(vx).(p.fr{ix});
		else
			val=arg;
		end
		if	numel(val) > 1
			sb='[]';
		end
			val={val};
% - cell
	case	'cell'
			sb='{}';
		if	par.hasvar
			val=arg(vx).(p.fr{ix});
		else
			val=arg;
		end

% - member
	otherwise
			par.level=par.level-1;
			return;

	end	% SWITCH

			[tpar,tpar]=run_vsize(1,p.fn{ix,3},p,par,...
				par.level,[p.fn{ix},sb],val{:});	%#ok
			par.ne=tpar.ne;
			par.b=tpar.b;
			par.r=tpar.r;
			par.level=par.level-1;

	end	% IF
end
%-------------------------------------------------------------------------------
function	par=show_ent(sflg,p,par,ix,varargin)
% general purpose print routine

		v=[p.fn(ix,:),varargin(1:4)];
		par.r=[par.r;v];
%D		disp(sprintf('SHOW> %5d %10s %s',sflg,varargin{[7,5]}));
	if	sflg
		s=sprintf(par.fmt,varargin{par.fmtb:end});
		disp(s);
	end
end
%-------------------------------------------------------------------------------