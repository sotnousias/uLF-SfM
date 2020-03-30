function [rayFeatures, rayDesc] = extractLfFeaturesInFrame(lf, minimumNumViews, matchThresh, peakThresh, onlyCentral, centId)
% extract light field SIFT-features in the same frame
% just for the 5x5 sub-aperture images

% imCent = lf{centralSubId, centralSubId};
% 
% [fCent, dCent] = vl_sift(single(rgb2gray(imCent)));

if (nargin < 3)
    matchThresh = 1.5; % default vlfeat value
end

if (nargin < 4)
    peakThresh = 0.0066;
end

if (nargin < 5)
    onlyCentral = 0;
end

if (nargin < 6)
    centId = (size(lf, 2) + 1) / 2;
end


lfFeatures = cell(size(lf));
lfDescriptors = cell(size(lf));



% get features and descriptos for all sub-aperture images
for x = 1:size(lf, 1)
    for y = 1:size(lf, 2)
        
%       more dense features
        [lfFeatures{x, y}, lfDescriptors{x, y}] = vl_covdet(im2single(rgb2gray(lf{x, y})), 'peakThreshold', peakThresh, 'MaxNumOrientations', 4);   
        lfDescriptors{x, y} = rootSift(lfDescriptors{x, y});
        
        
    end  
end

% central sub-aperture features and descriptors
fCent  = lfFeatures{centId, centId};
dCent = lfDescriptors{centId, centId};
numFeatures = size(fCent, 2);


if (onlyCentral)
    
    % for Illum
    centId = 3;
    centId = 3;
    rayDesc = cell(numFeatures, 1);
    rayFeatures = cell(numFeatures, 1);

    for mCent =  1 : size(fCent, 2)
        rayFeatures{mCent} = [fCent(1:2, mCent)' centId centId];
        rayDesc{mCent} = cat(1, rayDesc{mCent}, lfDescriptors{1, 1}(:, mCent)');
    end
    return;
end


rayFeatures = cell(numFeatures, 1);
descs = cell(numFeatures, 1);

for x = 1:size(lf, 1)
    for y = 1:size(lf, 2)
        
        if (x == centId && y == centId)
            continue;
        end
        
        matchIds = vl_ubcmatch(dCent, lfDescriptors{x, y}, matchThresh);
        
        for mId = 1:size(matchIds, 2)
           
            mCent = matchIds(1, mId);
            mQ = matchIds(2, mId);
            
            %check if already has been assigned a feature
            if isempty(rayFeatures{mCent})
                rayFeatures{mCent} = [fCent(1:2, mCent)' centId centId];
                descs{mCent} = dCent(:, mCent)';
            end
            
            rayFeatures{mCent} = cat(1, rayFeatures{mCent}, [lfFeatures{x, y}(1:2, mQ)' x y]);
            descs{mCent} = cat(1, descs{mCent}, lfDescriptors{x, y}(:, mQ)');
        end        
    end
end
   

rayFeaturesQ = rayFeatures;

rayFeatures = {};
rayDesc = {};

rayId = 1;

for n = 1:numel(rayFeaturesQ)

    if size(rayFeaturesQ{n}, 1) >= minimumNumViews
        rayFeatures{rayId} = rayFeaturesQ{n};
        rayDesc{rayId} = descs{n};
        rayId = rayId + 1;
    end
end





end

