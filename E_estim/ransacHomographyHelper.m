function [H, inl] = ransacHomographyHelper(matches1, matches2, K, th)
% Just convert pixels to homogeneous coordinates

if (nargin < 4)
    th = 1e-5;
end

if ~(size(matches1, 1)==2)
    matches1 = matches1';
end

if ~(size(matches2, 1)==2)
    matches2 = matches2';
end


x1 = K \ [matches1; ones(1, size(matches1, 2))];
x2 = K \ [matches2; ones(1, size(matches2, 2))];

[H, inl] = ransacfithomography(x1, x2, th);

end

