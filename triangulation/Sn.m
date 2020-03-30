% inlier selection with Rousseeuw & Croux's Sn scale estimator
% "Alternatives to the Median Absolute Deviation"
% Based on RousseeuwCrouxSn.m from http://www.ucl.ac.uk/~smgxprj/misc/RousseeuwCrouxSn.m

function inliers = Sn(X)
% Compute the measure of scale 'Sn', from Rousseeuw & Croux (1993)
%
%   A robust alternative to MADn for statistical outlier identification.
%   Unlike MADn, Sn does not make an assumption of symmetry, so in
%   principle should be more robust to skewed distributions.
%
%   The outputs of this function have been validated against equivalent
%   function in Maple(tm).
%
% Example:          X = [1 5 2 2 7 4 1 5]
%                   Sn = RousseeuwCrouxSn(X) % should give 3.015
%
%                   % use Sn to identify statistical outliers
%                	X = [1 5 2 2 7 50 1 5];
%                  	[Sn, x_j] = RousseeuwCrouxSn(X);
%                	outliers = X(x_j/Sn > 3) % criterion typically 2 or 3
%
% Requires:         none
%
% See also:         mad.m
%
% Author(s):        Pete R Jones <petejonze@gmail.com>
% 
% Version History:  19/04/2016	PJ  Initial version
%                                               
%
% Copyright 2016 : P R Jones
% *********************************************************************
% 

    % get number of elements
    n = length(X);

    % Set c: bias correction factor for finite sample size
    if n < 10
        cc = [NaN 0.743 1.851 0.954 1.351 0.993 1.198 1.005 1.131];
        c = cc(n);
    elseif mod(n,2)==0  % n is odd
        c = n/(n-.9);
    else                % n is even
        c = 1;
    end
    % compute median difference for each element
    x_j = nan(n,1);
    for i = 1:n
        x_j(i) = median(abs(X(i) - X([1:i-1 i+1:end]))); 
    end

    % compute median of median differences, and apply finite sample
    % correction, c, and 1.1926 correction for consistency with the
    % standard deviation over a Gaussian distribution
    Sn = 1.1926 * c * median(x_j);

    inliers = find(x_j <= 3.0*Sn);
end
