function [xim,yim,igrid,jgrid] = structvar2grid(structvar)
%function [xim,yim,igrid,jgrid] = structvar2grid(structvar)
%converts one of Ashely's famous structvars into a grid of points
%igrid,jgird <---> xim,yim

c1 = structvar(1:4:end,:);
c2 = structvar(2:4:end,:);
c3 = structvar(3:4:end,:);
c4 = structvar(4:4:end,:);

dy = mean([c4(:,2) - c1(:,2); c3(:,2) - c2(:,2)]);
dx = mean([c2(:,1) - c1(:,1); c3(:,1) - c4(:,1)]);

%min dist squared
mds = (0.25*dy).^2 + (0.25*dx).^2;

%iterate through and merge points close together to a single point
newsv = structvar;
keep = true([1 length(structvar)]);
for j = 1:length(structvar)
    if (keep(j))
        ds = structvar - repmat(structvar(j,:), [length(structvar), 1]);
        ds = sum(ds.^2, 2);
        close = find(ds < mds);
        
      %  plot (structvar(~keep,1), structvar(~keep,2), 'r.', structvar(keep,1), structvar(keep,2), 'b.', ...
       %     structvar(close,1), structvar(close,2), 'go', structvar(j,1), structvar(j,2), 'm*',...
        %    newsv(keep(1:j),1), newsv(keep(1:j), 2), 'y*');
        
 
       
        newsv(j,:) = mean(structvar(close,:),1);
        keep(close) = false;
        keep(j) = true;
    end
end
newsv = newsv(keep,:);
x0 = min(newsv(:,1));
y0 = min(newsv(:,2));

xim = newsv(:,1);
yim = newsv(:,2);
igrid = round((xim-x0)/dx);
jgrid = round((yim-y0)/dy);

%sort into rows and columns

[igrid,I] = sort(igrid);
jgrid = jgrid(I);
xim = xim(I);
yim = yim(I);

[jgrid,I] = sort(jgrid);
igrid = igrid(I);
xim = xim(I);
yim = yim(I);

        