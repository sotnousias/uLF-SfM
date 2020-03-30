% find inliers with the X84 rejection rule (MAD) for the outliers

function inliers = X84(res, n)
% n is the minimum number of inliers returned

  loc = median(res);
  res2 = abs(res-loc);
  scl = 1.4826*3.5*median(res2);
  inliers = find(res2<=scl);

  if length(inliers) < n % return the smallest n
    [y idx] = sort(res2);
    inliers = idx(1:n);
    %inliers = sort(inliers);
  end
end
