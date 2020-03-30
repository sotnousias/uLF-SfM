function [E, inliers] = ransacFitEssentialHelper(pts0, pts1, K0, K1, threshold)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

npts0=inv(K0)*[pts0; ones(1, size(pts0,2))];
npts1=inv(K1)*[pts1; ones(1, size(pts1,2))];

[E, inliers] = ransacfitessmatrix(npts0, npts1, threshold);


end

