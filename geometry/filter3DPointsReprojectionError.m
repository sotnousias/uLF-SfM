function [ids, reprojectionErrors] = filter3DPointsReprojectionError(pts3d, matchesAll, poses, subLfExt, K, reprTh)

scores = zeros(numel(pts3d), 1);

for n = 1:numel(pts3d)
    
    
   f1 = matchesAll(1).frame{n};
   f2 = matchesAll(2).frame{n};
   
   Rh = poses{f1}(1:3, 1:3);
   th = poses{f1}(1:3, 4);
   
   extMat1 = [Rh' -Rh'*th; zeros(1, 3) 1];
   Rh = poses{f2}(1:3, 1:3);
   th = poses{f2}(1:3, 4);
   extMat2 = [Rh' -Rh'*th; zeros(1, 3) 1];
     
   pt = pts3d(n).Location;
   
   ptH = [pt; 1];
   
   numViews1 = size(matchesAll(1).rays{n}, 1);
   numViews2 =  size(matchesAll(2).rays{n}, 1);
   
   reprs1 = zeros(numViews1, 1);
   reprs2 = zeros(numViews2, 1);
   
   for v = 1:numViews1
       
       subId = matchesAll(1).rays{n}(v, 3:4)';
       tSub = subLfExt{subId(1), subId(2)}(1:3, 4);
       tSub(2) = -tSub(2);
       impt =  matchesAll(1).rays{n}(v, 1:2)';
       repr = K * [eye(3) -tSub] * extMat1 * ptH;
       repr = repr ./ repr(3);
       repr = repr(1:2);       
       reprs1(v) = norm(impt-repr);
       
       
   end
   
   for v = 1:numViews2
       subId = matchesAll(2).rays{n}(v, 3:4)';
       tSub = subLfExt{subId(1), subId(2)}(1:3, 4);
       tSub(2) = -tSub(2);
       impt =  matchesAll(2).rays{n}(v, 1:2)';
       repr = K * [eye(3) -tSub] * extMat2 * ptH;
       repr = repr ./ repr(3);
       repr = repr(1:2);       
       reprs2(v) = norm(impt-repr);
       
       
   end
   
   sc1 = mean(reprs1);
   sc2 = mean(reprs2);
   
   scores(n) = (sc1 + sc2)/2; 
       
end

reprojectionErrors = scores;
   
ids = find(scores < reprTh);


end

