function [pts3d, pts] = lightFieldReconPtsRelativeFrames(matchesN, R1, t1, R2, t2, K, lf, lfId, getDesc, onlyCentral)
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
%
%
% lf                        -  Light field so that the color of the 3D
%                              point can be retrieved.
%
% lfId                      - light field id used to index the matches
%                            (for the colour of the point cloud).

if (nargin) < 9
    getDesc  = 0;
end

if (nargin < 10)
    onlyCentral = 0;
end
    
if (numel(matchesN(1).rays) ~= numel(matchesN(2).rays))
    error('Expecting same number of ray-bundle matches');
end

numPoints = numel(matchesN(1).rays);

% pts3d = zeros(3, numPoints);

extMat1 = [R1' -R1'*t1; zeros(1, 3) 1];
extMat2 = [R2' -R2'*t2; zeros(1, 3) 1];


for pt = 1:numPoints

    numMat1 = size(matchesN(1).rays{pt}, 1);
    numMat2 = size(matchesN(2).rays{pt}, 1);
    
    Ps = zeros(3, 4, numMat1 + numMat2);
    imPoints = zeros(2, numMat1 + numMat2);
    
    for n = 1:numMat1
        Ps(:, :, n) = K * [eye(3) [-matchesN(1).rays3d{pt}(n, 1:2)'; 0]] * extMat1;
        imPoints(:, n) = matchesN(1).rays{pt}(n, 1:2)';
    end
    
    for n =1:numMat2
        Ps(:, :, numMat1 + n) = K * [eye(3) [-matchesN(2).rays3d{pt}(n, 1:2)'; 0]] * extMat2;
        imPoints(:, numMat1 + n) = matchesN(2).rays{pt}(n, 1:2)';
    end

% %     debugging
%     fid = fopen(sprintf('point%d.txt', pt), 'w');
%     for n = 1:(numMat1+numMat2)
%         fprintf(fid, '%s %s\n', mat2str(squeeze(Ps(:, :, n))), mat2str(imPoints(:, n)));
%     end
%     fclose(fid);
%     
%     if (pt==79)
%         disp('Fack ...\n');
%     end
% pt
%     pt3d = nviewTriangulatePts(Ps, imPoints, [], 'dlt', 'vgg');
    [pt3d, inliers] = nviewTriangulatePts(Ps, imPoints, [383;552], 'rob', 'lm');
%     [pt3d, inliers] = nviewTriangulatePts(Ps, imPoints, [383;552], 'rob', 'lm');

%     [pt3d] = nviewTriangulatePts_v12(Ps, imPoints, [], 'rob', 'lm');


%     pt
    pts(:, pt) = pt3d;
    pts3d(pt).Location = pt3d;
    
    
    if isinf(pt3d)
        pts3d(pt).Color = inf(3, 1);
        pts3d(pt).inliers = [];
        continue;
    end
    
%     pts3d(pt).inliers = inliers;
    pts3d(pt).inliers = [];


    
    % keep the central descriptor for matching with sub-sequent frames
    if onlyCentral
        im = lf{1, 1};
    else
        im = lf{matchesN(lfId).rays{pt}(1, 3), matchesN(lfId).rays{pt}(1, 4)};
    end

    
pix = matchesN(lfId).rays{pt}(1, 1:2);
col = squeeze(im(round(pix(2)), round(pix(1)), :));
pts3d(pt).Color = col;
    
if (getDesc)
%     pts3d(pt).Descriptor = matchesN(1).descriptors{pt}(1, :);
%     pts3d(pt).Descriptor = extractFeatures(rgb2gray(im), pix, 'Method', 'SURF', 'FeatureSize', 128);
    % assign a random orientation for now
    fc = [pix'; 1; randn(1)];
    [feat, des] = vl_sift(im2single(rgb2gray(im)), 'frames', double(fc), 'Orientations');
    pts3d(pt).Descriptor = des(:, 1);
end
    
end


end

