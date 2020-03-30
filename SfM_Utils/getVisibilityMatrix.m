function [visMatrix, visBool] = getVisibilityMatrix(matches)

visMatrix = multimg_matches(matches);
visMatrix(visMatrix == -1) = 0;
visBool = visMatrix ~= 0;
end

