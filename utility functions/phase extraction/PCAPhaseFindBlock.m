function [thetaGlobal, coeffGlobal, reliabilityScoreGlobal, thetaInd, coeff, reliabilityScoreInd, ctr] = PCAPhaseFindBlock(datacube, nptsperramp)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

thetaInd = zeros(size(datacube,1), size(datacube,2), size(datacube,3)/nptsperramp);
reliabilityScoreInd = thetaInd;
%thetaGlobal = thetaInd;

%blocksize = 50;
chunksize = 1E8;
blocksize = floor(sqrt(chunksize / size(datacube,3)));
nrb = ceil(size(datacube, 1)/blocksize);
ncb = ceil(size(datacube, 2)/blocksize);

dr = (size(datacube,1) - blocksize)/(nrb - 1);
dc = (size(datacube,2) - blocksize)/(ncb - 1);

% ctr{1} = dc * (0:(ncb)-1) + blocksize/2;
% ctr{2} = dr * (0:(nrb)-1) + blocksize/2;

ctr = (0 + 0i) * zeros(1,nrb*ncb);
coeff = (0 + 0i) * zeros(nptsperramp, nrb*ncb);
ind = 1;
tic
for j = 1:nrb
    rinds = (1:blocksize) + floor((j-1)*dr);
    for k = 1:ncb
        ctr(ind) = dc * (k-1) + blocksize/2 + 1i*(dr * (j-1) + blocksize/2);
        cinds = (1:blocksize) + floor((k-1)*dc);
        rob = reshape(reshape(datacube(rinds,cinds,:), [], size(datacube,3)).',nptsperramp,[]).';
        [c,sc] = pca(rob, 'NumComponents', 2);
        magsq = sum(rob.^2,2);
        %z = sc(:,1) + 1i*sc(:,2);
        thetaInd(rinds,cinds,:) = permute(reshape(atan2(sc(:,2), sc(:,1)), [], blocksize, blocksize), [2 3 1]);
        reliabilityScoreInd(rinds,cinds,:) = permute(reshape(sum(sc.^2, 2)./magsq, [], blocksize, blocksize), [2 3 1]);

        coeff(:,ind) = c(:,1) + 1i*c(:,2);
        ind = ind+1;
    end
    disp ([num2str(j) ' rows of ' num2str(nrb) ' completed']);
    toc;
end
coeffGlobal = pca([real(coeff) imag(coeff)]', 'NumComponents', 2);
dcflat = reshape(reshape(datacube, [], size(datacube,3)).',20,[]).';
sc = dcflat * coeffGlobal; toc
magsq = sum(dcflat.^2,2);
reliabilityScoreGlobal = permute(reshape(sum(sc.^2, 2)./magsq, [], size(datacube,1), size(datacube,2)), [2 3 1]);
thetaGlobal = permute(reshape(atan2(sc(:,2), sc(:,1)), [], size(datacube,1), size(datacube,2)), [2 3 1]);
end

