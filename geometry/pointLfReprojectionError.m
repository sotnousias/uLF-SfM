function [reprs] = pointLfReprojectionError(pt, subMatches, pose, subLfExt, K)

ptH = [pt; 1];

Rh = pose(1:3, 1:3);
th = pose(1:3, 4);

extMat = [Rh' -Rh'*th; zeros(1, 3) 1];
numViews = size(subMatches, 1);

reprs = zeros(numViews, 1);

for v = 1:numViews
    
    subId = subMatches(v, 3:4)';
    tSub = subLfExt{subId(1), subId(2)}(1:3, 4);
    tSub(2) = -tSub(2);
    impt = subMatches(v, 1:2)';
    repr = K * [eye(3) -tSub] * extMat * ptH;
    repr = repr ./ repr(3);
    repr = repr(1:2);
    reprs(v) = norm(impt-repr);  
end

end
