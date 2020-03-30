function [posesBA, ptsCell] = bundleAdjustmentHelper(pts, visMatrix, poses, lfFeatures, registered, notRegistered, pts3dIds)
% writes the corresponding txt files to call BA.

[visMatrixOut, posesOut, lfFeaturesOut] = registeredPosesToSequential(visMatrix, poses, lfFeatures, registered, notRegistered);
% pointProjectionsToTxtV2(pts, pts3dIds, visMatrixOut, lfFeaturesOut, ptsStr);
% posesToTxt(posesOut, posesStr);
posesBA = posesToMat(posesOut);
ptsCell = pointProjectionsToCell(pts, pts3dIds, visMatrixOut, lfFeaturesOut);

end

