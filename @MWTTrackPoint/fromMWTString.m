function tp = fromMWTString (tp, str, camcalinfo)
%function tp = fromMWTString (tp, str, camcalinfo)
numSpinePts = 11;


% expr = [numtok('frame') numtok('et') numtok('comx') numtok('comy') numtok('area') numtok('v1x') numtok('v1y') numtok('cov2') numtok('len') numtok('wid') ...
%     '(%\s*)(?<spinestr>[-,\d,\.,\s]*)(%%\s*)' numtok('ctrx0') numtok('ctry0') numtok('ncpts') '(?<ctrstr>[\S]*)']
expr = '(?<frame>[-,\d,\.]+)(\s+)(?<et>[-,\d,\.]+)(\s+)(?<comx>[-,\d,\.]+)(\s+)(?<comy>[-,\d,\.]+)(\s+)(?<area>[-,\d,\.]+)(\s+)(?<v1x>[-,\d,\.]+)(\s+)(?<v1y>[-,\d,\.]+)(\s+)(?<cov2>[-,\d,\.]+)(\s+)(?<len>[-,\d,\.]+)(\s+)(?<wid>[-,\d,\.]+)(\s+)(%\s*)(?<spinestr>[-,\d,\.,\s]*)(%%\s*)(?<ctrx0>[-,\d,\.]+)(\s+)(?<ctry0>[-,\d,\.]+)(\s+)(?<ncpts>[-,\d,\.]+)(\s+)(?<ctrstr>[\S]*)';

bob = regexp(str, expr, 'names');

if (isempty(bob))
    [nums,~,~,ind] = sscanf(str, '%g', 10);

    tp.ind = nums(1);
    tp.et = nums(2);

    tp.loc = [nums(3);nums(4)];

    tp.area = nums(5);
    % 
    % 
    % v1 = [num(6);num(7)];
    % v2 = [num(8);num(9)];
    % d1 = sqrt(sum(v1.^2));
    % d2 = sqrt(sum(v2.^2));
    % 
    % v1 = v1/d1;
    % v2 = v2/d2;
    ind = ind + find(str(ind:end) == '%', 1,'first');
    [spn,~,~,ind2] = sscanf(str(ind:end), '%g', [2 numSpinePts]);
    spn(1,:) = spn(1,:) + tp.loc(1);
    spn(2,:) = spn(2,:) + tp.loc(2);
    tp.spine = spn;

    ind = ind + ind2;
    ind = ind + find(str(ind:end) == '%', 1,'first') + 1;

    [ctr0,~,~,ind2] = sscanf(str(ind:end), '%g', [2 1]);
    ind = ind + ind2;
    [npts,~,~,ind2] = sscanf(str(ind:end), '%d', 1);
    ind = ind + ind2;
    ccstr = regexp(str(ind:end), '\S*', 'match', 'once');

    [x,y] = unpackKerrChainCode(ccstr, npts);
    x = x+ctr0(1);
    y = y+ctr0(2);
    tp.contour = [x;y];
else
    tp.ind = sscanf(bob.frame,'%d');
    tp.loc = [sscanf(bob.comx,'%g');sscanf(bob.comy,'%g')];
    tp.et = sscanf(bob.et,'%g');
    tp.area = sscanf(bob.area,'%g');
    spn = reshape(sscanf(bob.spinestr, '%g'),2,[]);
    spn(1,:) = spn(1,:) + tp.loc(1);
    spn(2,:) = spn(2,:) + tp.loc(2);
    tp.spine = spn;
    ctr0 = [sscanf(bob.ctrx0,'%g');sscanf(bob.ctry0,'%g')];
    npts = sscanf(bob.ncpts,'%d');
    [x,y] = unpackKerrChainCode(bob.ctrstr, npts);
    x = x+ctr0(1);
    y = y+ctr0(2);
    tp.contour = [x;y];
    v1 = [sscanf(bob.v1x,'%g');sscanf(bob.v1y,'%g')];
    cov1 = sqrt(sum(v1.^2));
    v = v1./cov1;
    D = [cov1,0;0,sscanf(bob.cov2,'%g')];
    V = [v,[-v(2);v(1)]];
    c = V*D/V;
    tp.cov = [c(1) c(2) c(4)];
end

if (~isempty(camcalinfo))
     tp.loc = cameracalibration.realPtsFromCamPts(tp.loc);
     tp.spine = cameracalibration.realPtsFromCamPts(tp.spine);
     tp.contour = cameracalibration.realPtsFromCamPts(tp.contour);
end

tp.head = tp.spine(:,end);
tp.tail = tp.spine(:,1);
tp.mid = tp.spine(:,ceil(size(tp.spine,2)/2));
tp.htValid = true;

function patt = numtok(tokenName)
patt = ['(?<' tokenName '>[-,\d,\.]+)(\s+)'];
