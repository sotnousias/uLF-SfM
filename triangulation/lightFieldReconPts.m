function [pts3d] = lightFieldReconPts(matchesN, R, t, K)
% matchesN               -   Light field matches struct with the following fields:
% 
% matchesN(i).rays{n}       -   Nx4 matrix where the first two columns are the pixel and 
%                               the last two the index of the sub-aperture image the pixel 
%                               belongs. The index i corresponds to the
%                               light field frame.
% matchesN(i).rays3d{n}     -  Nx4 matrix where the first two columns are
%                              the camera center (z-coordinate is zero) and the 
%                              last two the homogeneous x-y coordinates of a 
%                              pixel (i.e. multiplied by the inverse of K). 



if (numel(matchesN(1).rays) ~= numel(matchesN(2).rays))
    error('Expecting same number of ray-bundle matches');
end

numPoints = numel(matchesN(1).rays);

pts3d = zeros(3, numPoints);

extMat = [R' -R'*t; zeros(1, 3) 1];


for pt = 1:numPoints

    numMat1 = size(matchesN(1).rays{pt}, 1);
    numMat2 = size(matchesN(2).rays{pt}, 1);
    
    Ps = zeros(3, 4, numMat1 + numMat2);
    imPoints = zeros(2, numMat1 + numMat2);
    
    for n = 1:numMat1
        Ps(:, :, n) = K * [eye(3) [-matchesN(1).rays3d{pt}(n, 1:2)'; 0]];
        imPoints(:, n) = matchesN(1).rays{pt}(1:2)';
    end
    
    for n =1:numMat2
        Ps(:, :, numMat1 + n) = K * [eye(3) [-matchesN(2).rays3d{pt}(n, 1:2)'; 0]] * extMat;
        imPoints(:, numMat1 + n) = matchesN(2).rays{pt}(1:2)';
    end

    pts3d(:, pt) = nviewTriangulatePts(Ps, imPoints, [], 'dlt', 'lm');
    
    
end


end

