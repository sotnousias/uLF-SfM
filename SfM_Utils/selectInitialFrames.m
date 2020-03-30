function [f1, f2] = selectInitialFrames(geomVer, numMatches)
%
geomVerPairMatch = geomVer .* numMatches;
[val, id] = max(geomVerPairMatch(:));

[f2, f1] = ind2sub(size(geomVer), id);

end

