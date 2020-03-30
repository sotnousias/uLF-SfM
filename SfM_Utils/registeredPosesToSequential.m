function [visMatrixOut, posesOut, lfFeaturesOut] = registeredPosesToSequential(visMatrix, poses, lfFeatures, registered, notRegistered)
% Map the registered poses to that of a sequential index of the cameras

n = numel(registered);

lfFeaturesOut = cell(n ,1);
posesOut = cell(n, 1);

% we dont want to change the 1st dimension as we are indexing the
% reconstructed 3D points.
visMatrixOut = zeros(size(visMatrix, 1), n);

for n = 1:numel(registered)
    posesOut{n} = poses{registered(n)};
    lfFeaturesOut{n} = lfFeatures{registered(n)};
    visMatrixOut(:, n) = visMatrix(:, registered(n));
end

end

