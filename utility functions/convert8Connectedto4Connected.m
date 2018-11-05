function [dx4,dy4] = convert8Connectedto4Connected(dx8, dy8)
% function [dx4,dy4] = convert8Connectedto4Connected(dx8, dy8)
% converts and 8-connected chain to a 4 connected chain

dxmat = {[-1 0], [0], [1 0]; [-1], [], [1]; [-1 0], [0], [1 0]};
dymat = {[0 -1], [-1], [0 -1]; [0], [], [0]; [0 1], [1], [0 1]};

dx4 = [dxmat{sub2ind(size(dxmat), dy8+2,dx8+2)}];
dy4 = [dymat{sub2ind(size(dymat), dy8+2,dx8+2)}];
