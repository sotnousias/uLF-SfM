function [pts3d, pts3D, ids, reprojectionErrs, pts3dIds, pts, cols, pc] = initialiseStructure(matchesN, poses, f1, f2, K, subLfExt, lf, reprTh, pts3dIds, lfId, getDesc)

if (nargin < 10)
    lfId = 2;
end

if (nargin < 11)
    getDesc = 0;
end
%

R1 = poses{f1}(1:3, 1:3);
t1 = poses{f1}(1:3, 4);
R2 = poses{f2}(1:3, 1:3);
t2 = poses{f2}(1:3, 4);

[pts3d, pts3D] = lightFieldReconPtsRelativeFrames(matchesN, R1, t1, R2, t2, K, lf, lfId, getDesc);

matchesAll = matchesN;

for n = 1:numel(matchesN(1).rays)
    
    matchesAll(1).frame{n} = f1;
    matchesAll(2).frame{n} = f2;
end

% posesN{1} = poses{f1};
% posesN{2} = poses{f2};

% if you filter the outliers we need to remove the ptIds as well

%
[ids, reprojectionErrs] = filter3DPointsReprojectionError(pts3d, matchesAll, poses, subLfExt, K, reprTh);


pts = [];
cols = [];
for n = 1:numel(pts3d)
    
    pts = [pts pts3d(n).Location];
    cols = [cols; pts3d(n).Color'];
end

pts = pts(:, ids);
cols = cols(ids, :);
pts3dIds = pts3dIds(ids);
pc = pointCloud( pts', 'Color', uint8(cols));

end

