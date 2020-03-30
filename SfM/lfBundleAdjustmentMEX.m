function [poses, pts] = lfBundleAdjustmentMEX(intrStr, lfIntrStr, ptsCell, posesMat, poses, registered)
% We assume that the bundleAdjustmentHelper has already been called
% already.

% next two lines should be moved elsewhere and called only once!
K = textread(intrStr);
subRelPosesMat = textread(lfIntrStr);
[posesTmp, pts] = eucsbademo(posesMat, subRelPosesMat, ptsCell, K); % MEX

pts = pts';

for n = 1:numel(registered)
    Rh = quat2rot(posesTmp(n, 1:4));
    th = posesTmp(n, 5:end)';
    poses{registered(n)} = [Rh th];
end

end

