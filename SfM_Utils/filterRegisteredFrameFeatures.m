function [visMatrix, visBool, lfFeatures] = filterRegisteredFrameFeatures(visMatrix, visBool, pts, pts3dIds, poses, frameId, lfFeatures, K, subLfExt)
% filter features of already reconstructed points based on reprojection
% error

outlReprTh = 1.5;

ptIds = find(visMatrix(pts3dIds, frameId));
pts = pts(:, ptIds);
pts3dIds = pts3dIds(ptIds);

for n = 1:numel(pts3dIds)
    
   feats = lfFeatures{frameId}{visMatrix(pts3dIds(n), frameId)};
   errs = pointLfReprojectionError(pts(:, n), feats, poses{frameId}, subLfExt, K) ;
   
   rmIds = find(errs > outlReprTh);
   
   if isempty(rmIds)
       continue;
   end
   
   if (numel(rmIds) == size(feats, 1))
       visMatrix(pts3dIds(n), frameId) = 0;
       visBool(pts3dIds(n), frameId) = 0;
   else
       
       feats(rmIds, 1) = inf(numel(rmIds), 1);
       feats = feats(~(feats(:, 1)==inf), :);
       lfFeatures{frameId}{visMatrix(pts3dIds(n), frameId)} = feats;
   end
    
end

