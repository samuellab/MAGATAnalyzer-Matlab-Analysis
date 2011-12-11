function [x,y] = unpackKerrChainCode(ccstr, npts)

dx = zeros([1, npts]);
dy = dx;
ccstr = uint8(ccstr) - '0'; %0 = 48 not 33, as described in the manual
bc = [ccstr ccstr ccstr];
bc(1:3:end) = bitand(bitshift(ccstr, -4), 3);
bc(2:3:end) = bitand(bitshift(ccstr, -2), 3);
bc(3:3:end) = bitand(ccstr, 3);

bc = bc(1:(npts));
dx(bc == 0) = -1;
dx(bc == 1) = 1;
dy(bc == 2) = -1;
dy(bc == 3) = 1;

x = [0 cumsum(dx)];
y = [0 cumsum(dy)];
% 
% 
% for j = 2:npts
%     c = ccstr(floor((j+1)/3)); 
%     bs = 4 - 2*mod((j+1),3);
%     %cc(j) = c;
%     x(j) = x(j-1);
%     y(j) = y(j-1);
%     switch (bitand(bitshift(c, -bs), 3))
%         case 0
%             x(j) = x(j) - 1;
%         case 1
%             x(j) = x(j) + 1;
%         case 2
%             y(j) = y(j) - 1;
%         case 3
%             y(j) = y(j) + 1;
%         otherwise
%             disp ('you dumb!');
%     end
%     %plot (x(1:j), y(1:j), 'b.-'); pause
% end
% ccstr
% char(cc + 48)