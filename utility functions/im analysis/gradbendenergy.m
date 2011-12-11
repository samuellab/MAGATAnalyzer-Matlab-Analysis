function gb = gradbendenergy(ctr,closed)
existsAndDefault('closed', true);
if (closed)
    dx = 0.5 * (diff(ctr(:, [end 1:end]),[],2) + diff(ctr(:, [1:end 1]),[],2));
    d4x = diff(ctr(:, [end-1:end 1:end 1:2]), 4, 2);
else    
   % ctr = ctr(:, [1 1 1:end end end]);
   % ctr = [2*ctr(:,[1 1])-ctr(:,[3 2]), ctr, 2*ctr(:, [end end])-ctr(:, [end-1 end-2])];   
    firstpoint = ctr(:,3) + 4*(ctr(:,3)-ctr(:,4));
    secondpoint = 1/2 * (firstpoint + ctr(:,1));
    lastpoint = ctr(:,end-2) + 4*(ctr(:,end-2) - ctr(:,end-3));
    secondtolastpoint = 1/2*(lastpoint + ctr(:,end));
%    ctr = [interp1(ctr', [-1 0], 'spline', 'extrap')', ctr, interp1(ctr', [length(ctr) + (1:2)], 'spline', 'extrap')'];
    newctr = [firstpoint secondpoint ctr secondtolastpoint lastpoint];
    dx = 0.5 * (diff(newctr(:, 2:(end-2)),[],2) + diff(newctr(:, 3:(end-1)),[],2));
    d4x = diff(newctr, 4, 2);   
    
    %make forces at the end perpendicular to tangent
    that = diff(ctr(:,[1 2]),[],2); that = that./sqrt(sum(that.^2));
   % whos that
   % size(sum(that.*d4x(:,1)))
    d4x(:,1) = d4x(:,1) - that*sum(that.*d4x(:,1));
    that = diff(ctr(:,[2 3]),[],2); that = that./sqrt(sum(that.^2));
    d4x(:,2) = d4x(:,2) - that*sum(that.*d4x(:,2));
    that = diff(ctr(:,[end-1 end]),[],2); that = that./sqrt(sum(that.^2));
    d4x(:,end) = d4x(:,end) - that*sum(that.*d4x(:,end));
    that = diff(ctr(:,[end-2 end-1]),[],2); that = that./sqrt(sum(that.^2));
    d4x(:,end-1) = d4x(:,end-1) - that*sum(that.*d4x(:,end-1));

end
dl = sqrt(sum(dx.^2));
denom = [sum(dx.^2);sum(dx.^2)].^(1.5);
ds = dl/sum(dl) + eps;
%plot (ds); hold all
%plot (d4x')
%gb = d4x./([ds;ds].^2);
gb = d4x./denom.*[ds;ds];
