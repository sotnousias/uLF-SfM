function [Rw, tw, registerSucc, ptsW, rays3, idsq, inlAbs] = registerLfAbsolutePose(visBool, visMatrix, pts3dIds, nextFrame, pts, lfFeat, K, subLfExt, ransacAbsTh, maxIters)
%

if (nargin < 9)
    ransacAbsTh = 0.25;
end

if (nargin < 10)
    maxIters = 20;
end

seenIds = find(visBool(pts3dIds, nextFrame));

Rw = zeros(3);
tw = zeros(3, 1);
registerSucc = 0;

ptsW = [];
rayFeats3 = {};
idsq = [];

for n = 1 : numel(seenIds)
    
    % which point
    ptId = find( pts3dIds == pts3dIds(seenIds(n)) );
    % which feature of the to-be-registered frame
    featId = visMatrix(pts3dIds(ptId), nextFrame);
    %
    
    %     ptId = pts3dIds(seenIds(n));
    %     featId = visMatrix(pts3dIds(seenIds(n)), nextFrame);
    
    %     pt = pts3d(ptId).Location;
    pt = pts(:, ptId);
    
    numRays = size(lfFeat{featId}, 1);
    rayFeats3.rays3d{n} = lfPixelsToRays(lfFeat{featId}, K, subLfExt) ;
    ptsW = [ptsW repmat(pt, [1 numRays])];
    
    % the follwing is [pt3dId numFeatures ptId] where pt3dId is the index in
    % the visMatrix and ptId is the index in the pts array
    % this makes easier to remove outliers since it is the same size like
    % rays3
    tmp = [repmat(pts3dIds(ptId), [numRays 1]) (1:numRays)' repmat(ptId, [numRays 1])];
    idsq = [idsq; tmp];
    
end

rays3 = rayCellToMatrix(rayFeats3.rays3d);
rays3 = prepareRaysForOpenGV(rays3');


[X, inlAbs] = opengvV2('gp3p_ransac', ptsW, rays3');
%     [X, inlAbs] = opengv('gp3p_ransac', ptsW, rays3');

fprintf('frame %d: %.2f inliers\n', nextFrame, numel(inlAbs)/size(rays3, 1));

iter = 1;

while (numel(inlAbs)/size(rays3, 1) < ransacAbsTh) && iter<maxIters
    
    [X, inlAbs] = opengvV2('gp3p_ransac', ptsW, rays3');
%     fprintf('frame %d: %.2f inliers\n', nextFrame, numel(inlAbs)/size(rays3, 1));
    iter = iter + 1;
    
end

if (iter==maxIters && (numel(inlAbs)/size(rays3, 1) < ransacAbsTh))
    return;    
end

Xout = opengv('abs_nonlin_noncentral', double(inlAbs), ptsW, rays3', X);

Rw = Xout(1:3, 1:3);
tw = Xout(1:3, 4);
registerSucc = 1;



end

