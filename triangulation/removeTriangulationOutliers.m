function [lfFeatures, matchesAll] = removeTriangulationOutliers(lfFeatures, matchesAll, pts3d, regId, nextFrame, newIds, visMatrix)
% The assumption here is that regId is frame 1 and nextFrame is frame 2

for n = 1:numel(newIds)
    
    if isinf(pts3d(n).Location)
        continue;
    end
    
    if isempty(pts3d(n).inliers)
        continue;
    end
    inliers = pts3d(n).inliers;
    numMat1 = size(lfFeatures{regId}{visMatrix(newIds(n), regId)}, 1);
    tmp = lfFeatures{regId}{visMatrix(newIds(n), regId)};
        
    ids = find(inliers <= numMat1);
    ids2 = find(inliers > numMat1);
    
    if (isempty(ids) || isempty(ids2))
        continue;
    end
    
    tmp = tmp(inliers(ids), :);
    
    lfFeatures{regId}{visMatrix(newIds(n), regId)} = tmp;
    tmp = matchesAll(1).rays{n};
    tmp = tmp(inliers(ids), :);
    matchesAll(1).rays{n} = tmp;
    tmp = matchesAll(1).rays3d{n};
    tmp = tmp(inliers(ids), :);
    matchesAll(1).rays3d{n} = tmp;
    
    ids = find(inliers > numMat1);
    tmp = lfFeatures{nextFrame}{visMatrix(newIds(n), nextFrame)};
    tmp = tmp(inliers(ids)-numMat1, :);
    
    lfFeatures{nextFrame}{visMatrix(newIds(n), nextFrame)} = tmp;
    tmp = matchesAll(2).rays{n};
    tmp = tmp(inliers(ids)-numMat1, :);
    matchesAll(2).rays{n} = tmp;
    tmp = matchesAll(2).rays3d{n};
    tmp = tmp(inliers(ids)-numMat1, :);
    matchesAll(2).rays3d{n} = tmp;
end
    
    


end

