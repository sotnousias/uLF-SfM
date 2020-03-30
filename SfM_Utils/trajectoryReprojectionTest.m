function [isOutlier, ids, reprojectionErrors] = trajectoryReprojectionTest(pt3d, feats, framePose, subLfExt, K, reprTh)

reprs = zeros(size(feats, 1), 1);

Rh = framePose(1:3, 1:3);
th = framePose(1:3, 4);

extMat1 = [Rh' -Rh'*th; zeros(1, 3) 1];

for n = 1:numel(reprs)

%    pt = pts3d(n).Location;
   
   ptH = [pt3d; 1];
   
   subId = feats(n, 3:4)';
   tSub = subLfExt{subId(1), subId(2)}(1:3, 4);
   tSub(2) = -tSub(2);
   impt =  feats(n, 1:2)';
   repr = K * [eye(3) -tSub] * extMat1 * ptH;
   repr = repr ./ repr(3);
   repr = repr(1:2);
   reprs(n) = norm(impt-repr);
       
end

reprojectionErrors = reprs; 
ids = find(reprs < reprTh);
outl = find(reprs >= reprTh);

if numel(outl) == numel(reprs)
    isOutlier = true;
else
    isOutlier = false;
    
end


end

