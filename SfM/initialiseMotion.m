function [visMatrix, visBool, lfFeatures, matchesNfil, pts3dIds, R, t] = initialiseMotion(lfFeatures, visBool, visMatrix, K, subLfExt, centDes, f1, f2, ransacRelTh, sameSubviews)

if (nargin < 9)
    ransacRelTh = 0.8;
%     ransacRelTh = 0.35;

end
if (nargin < 10)
    sameSubviews = 0;
end

ids = visBool(:, f1) .* visBool(:, f2);
ids = find(ids);
pts3dIds = ids;

for n = 1:numel(ids)

matchesN(1).rays{n} = lfFeatures{f1}{visMatrix(ids(n), f1)};
matchesN(1).descriptors{n} = centDes{f1}(:, visMatrix(ids(n), f1));

matchesN(2).rays{n} = lfFeatures{f2}{visMatrix(ids(n), f2)};
matchesN(2).descriptors{n} = centDes{f2}(:, visMatrix(ids(n), f2));

end

for n = 1:numel(matchesN(1).rays)
    
    matchesN(1).rays3d{n} = lfPixelsToRays(matchesN(1).rays{n}, K, subLfExt);
    
    matchesN(2).rays3d{n} = lfPixelsToRays(matchesN(2).rays{n}, K, subLfExt);
        
end

[rays1, rays2, ids1, ids2] = getRayMatrices(matchesN, sameSubviews);
%
rays1 = prepareRaysForOpenGV(rays1);
rays2 = prepareRaysForOpenGV(rays2);


%
% percentage of ransac inliers
ransacRelTh = 0.75;

[X, inliers17pt] = opengvV2('seventeenpt_ransac', rays1', rays2');
fprintf('Inliers: %f...\n', numel(inliers17pt) / size(rays1, 1));


while (numel(inliers17pt)/size(rays1, 1) < ransacRelTh)
    [X, inliers17pt] = opengvV2('seventeenpt_ransac', rays1', rays2');
%     X
    fprintf('Inliers: %f...\n', numel(inliers17pt) / size(rays1, 1));
%     pause();

end

[visMatrix, visBool, lfFeatures, matchesNfil, outlierPtIds] = removeOutliersLfSfMRelPose(matchesN, inliers17pt, ids1, ids2, f1, f2, pts3dIds, lfFeatures, visMatrix, visBool);

% looks like there is an index mismatch bewtween ids
pts3dIds = pts3dIds(~ismember(1:numel(pts3dIds), outlierPtIds));

%
Xout = opengv('rel_nonlin_noncentral', double(inliers17pt), rays1', rays2', X);
Rgv = Xout(1:3, 1:3);
tgv = Xout(1:3, 4);
% not use if I should use this for now?
% tgv = Xout(1:3, 4) ./ norm(Xout(1:3, 4)) * ( norm(X(1:3, 4)) + norm(Xout(1:3, 4)))/2;

%
R = Rgv;
t = tgv;


end

