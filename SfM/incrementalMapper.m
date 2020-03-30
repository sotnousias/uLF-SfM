function [visBool, visMatrix, poses, pts, cols, pts3dIds, lfFeatures] = incrementalMapper(lf, visBool, visMatrix, regId, nextFrame, lfFeatures, centDes, K, subLfExt, poses, pts3dIds, pts, cols, reprTh, registered )

matchIds = visBool(:, regId) .* visBool(:, nextFrame);
matchIds = find(matchIds);

newIds = matchIds(~ismember(matchIds, pts3dIds));

if (isempty(newIds))
    return;
end

for n = 1:numel(newIds)
    ptId = newIds(n);
    matchesN(1).rays{n} = lfFeatures{regId}{visMatrix(ptId, regId)};
    matchesN(1).descriptors{n} = centDes{regId}(:, visMatrix(ptId, regId));
    matchesN(2).rays{n} = lfFeatures{nextFrame}{visMatrix(ptId, nextFrame)};
    matchesN(2).descriptors{n} = centDes{nextFrame}(:, visMatrix(ptId, nextFrame));
    
end

for n = 1:numel(matchesN(1).rays)
    matchesN(1).rays3d{n} = lfPixelsToRays(matchesN(1).rays{n}, K, subLfExt);
    matchesN(2).rays3d{n} = lfPixelsToRays(matchesN(2).rays{n}, K, subLfExt);
end

R1 = poses{regId}(1:3, 1:3);
t1 = poses{regId}(1:3, 4);
R2 = poses{nextFrame}(1:3, 1:3);
t2 = poses{nextFrame}(1:3, 4);

[pts3d, pts3D] = lightFieldReconPtsRelativeFrames(matchesN, R1, t1, R2, t2, K, lf, 2);

npts = [];
ncols = [];
for n = 1:numel(pts3d)
    npts = [npts pts3d(n).Location];
    ncols = [ncols; pts3d(n).Color'];
end

% filter points based on reprojection
matchesAll = matchesN;

for n = 1:numel(matchesN(1).rays)
    
    matchesAll(1).frame{n} = regId;
    matchesAll(2).frame{n} = nextFrame;
end

[ids, reprojectionErrs] = filter3DPointsReprojectionError(pts3d, matchesAll, poses, subLfExt, K, reprTh);

if isempty(ids)
    return;
end

for n = 1:numel(ids)
    
    % id to index the visMatrix
    pt3dId = newIds(ids(n));
    
    % we need only already registered frames and not current ones
    % as well
    qFrames = find(visBool(pt3dId, :));
    qFrames = qFrames(ismember(qFrames, registered));
    qFrames = qFrames(~ismember(qFrames, [regId nextFrame]));
    
    if ~isempty(qFrames)
        
        for k = qFrames
            
            [isOutlier, inlIds, reprs] = trajectoryReprojectionTest(npts(:, ids(n)),...
                lfFeatures{k}{visMatrix(pt3dId, k)}, poses{k}, subLfExt, K, reprTh);
            
            if isOutlier
                visMatrix(pt3dId, k) = 0;
            else
                % keep only the inliers
                lfFeatures{k}{visMatrix(pt3dId, k)} = lfFeatures{k}{visMatrix(pt3dId, k)}(inlIds, :);
                
            end
        end
    end
end

pts3dIds = [pts3dIds; newIds(ids)];

if (find(pts3dIds == 4202))
    5;
end

pts = [pts npts(:, ids)];
cols = [cols;  ncols(ids, :)];
        
end

