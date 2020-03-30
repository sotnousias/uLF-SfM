function [visMatrix, visBool, lfFeatures, matchesFiltered, outlierPtIds] = removeOutliersLfSfMRelPose(matchesN, inliers, ids1, ids2, f1, f2, pts3dIds, lfFeatures, visMatrix, visBool, numOutlViews, debugFlag)
% Remove the outliers from the matchesN struct. Returned are the
% outliers-removed matches as well as the ids of the 3D points removed.
%%%%%%%%%% TO DO: provide lf as input for debugging purposes %%%%%%

if (nargin < 11)
    % if one feature is outlier, it will have all the rays (assumed
    % correct) in the other frame as mismatches, we relax this since we
    % may have in frame mis-matches in both light fields to be registered. 
    numOutlViews = 4;
end


if (nargin < 12)
    debugFlag = 0;
end

% [nx, ny, nz] = size(lf{1}{3, 3});

n = 1:size(ids1, 1);
ids = n(~ismember(n, inliers));
outl = ids;

% ids of the outlier ray-pairs
outIds1 = ids1(outl, :);
outIds2 = ids2(outl, :);

curId = 0;
ptOutliers = 0;
featOutliers = 0;

outlierPtIds = [];

while (curId < size(outIds1, 1))

    ptId = outIds1(curId + 1, 1);
    % this is the same for the outIds2 as well
    qids = find(outIds1(:, 1) == ptId);
    
    % check if this point is a total outlier, i.e. discard all the matches
    if numel(qids) == (size(matchesN(1).rays{ptId}, 1) * size(matchesN(2).rays{ptId}, 1))
        fprintf('Point %d is an outlier...\n', ptId);
        ptOutliers = ptOutliers + 1;
        matchesN(1).rays{ptId} = {};
        matchesN(2).rays{ptId} = {};
        matchesN(2).rays3d{ptId} = {};
        matchesN(2).rays3d{ptId} = {};
        outlierPtIds = [outlierPtIds ptId];
        
    else % we may still have some mis-detected features in sub-aperture images
        numRaysL = size(matchesN(1).rays{ptId}, 1);
        numRaysR = size(matchesN(2).rays{ptId}, 1);
        
        % now find the number a specific feature is observed and if it macthes
        % the number of rays on the other lf, remove it, otherwise do nothing,
        % i.e. it was selected as an outlier due to strict threshold
        featIds = unique(outIds1(qids, 2));
        for n = 1:numel(featIds)
            cs = find(outIds1(qids, 2) == featIds(n));
            if numel(cs) > numOutlViews % == numRaysR
                if (debugFlag)
                    fprintf('Found left frame outlier at point %d, id: %d\n', ptId, featIds(n));
                end
                featOutliers = featOutliers + 1;
%                 if (debugFlag)
%                     visualiseLfFeatures(lf{5}, matchesN(1).rays{pt});
%                     hold on;
%                     pix(1) = matchesN(1).rays{pt}(featIds(n), 1) + (matchesN(1).rays{pt}(featIds(n), 4) - 1)*ny;
%                     pix(2) = matchesN(1).rays{pt}(featIds(n), 2) +  (matchesN(1).rays{pt}(featIds(n), 3) - 1)*nx;
%                     plot(pix(1), pix(2), 'og')
%                     pause();
%                     close all;
%                 end
                matchesN(1).rays{ptId}(featIds(n), :) = inf(1, 4);
            end
        end
        
        % now check the right light-field
        featIds = unique(outIds2(qids, 2));
        for n = 1:numel(featIds)
            cs = find(outIds2(qids, 2) == featIds(n));
            if numel(cs) > numOutlViews %== numRaysL
                if (debugFlag)
                    fprintf('Found right frame outlier at point %d, id: %d\n', ptId, featIds(n));
                end
                featOutliers = featOutliers + 1;
%                 if (debugFlag)
%                     visualiseLfFeatures(lf{2}, matchesN(2).rays{pt});
%                     hold on;
%                     pix(1) = matchesN(2).rays{pt}(featIds(n), 1) + (matchesN(2).rays{pt}(featIds(n), 4) - 1)*ny;
%                     pix(2) = matchesN(2).rays{pt}(featIds(n), 2) +  (matchesN(2).rays{pt}(featIds(n), 3) - 1)*nx;
%                     plot(pix(1), pix(2), 'og')
%                     pause();
%                     close all;
%                 end
                matchesN(2).rays{ptId}(featIds(n), :) = inf(1, 4); 
            end
        end
    end
    curId = curId + numel(qids);
end

clear id;
id = 1;
for n = 1:numel(matchesN(1).rays)
%     n
    if ( isempty(matchesN(1).rays{n}) || isempty(matchesN(2).rays{n}) )
        continue;
    else
        % left frame
        f1s = matchesN(1).rays{n};
        goodIds1 = find(~(f1s(:, 1) == inf));
        % right frame
        f2s = matchesN(2).rays{n};
        goodIds2 = find(~(f2s(:, 1) == inf));
        % ensure no ray-feature is empty in both frames so that we can have
        % valid triangulation 
        if (isempty(goodIds1) || isempty(goodIds2))
            % update outlier points counter and index
            outlierPtIds = [outlierPtIds n];
            ptOutliers = ptOutliers + 1;
            continue;
        end
        
        matchesOut(1).rays{id} = matchesN(1).rays{n}(goodIds1, :);
        matchesOut(1).rays3d{id} = matchesN(1).rays3d{n}(goodIds1, :);
        matchesOut(1).descriptors{id} = matchesN(1).descriptors{n};
        
        matchesOut(2).rays{id} = matchesN(2).rays{n}(goodIds2, :);
        matchesOut(2).rays3d{id} = matchesN(2).rays3d{n}(goodIds2, :);
        matchesOut(2).descriptors{id} = matchesN(2).descriptors{n};
        
        % now clear the lfFeatures as well, for Bundle Adjustment
        ptId = pts3dIds(n);
        lfFeatures{f1}{visMatrix(ptId, f1)} = matchesOut(1).rays{id};
        lfFeatures{f2}{visMatrix(ptId, f2)} = matchesOut(2).rays{id};
        
    end
    id = id + 1;
end

fprintf('Outlier points: %d, outlier features: %d\n', ptOutliers, featOutliers);

matchesFiltered = matchesOut;
outlierPtIds = sort(outlierPtIds);

visMatrix(pts3dIds(outlierPtIds), f1) = zeros(numel(outlierPtIds), 1);
visBool(pts3dIds(outlierPtIds), f1) = zeros(numel(outlierPtIds), 1);

visMatrix(pts3dIds(outlierPtIds), f2) = zeros(numel(outlierPtIds), 1);
visBool(pts3dIds(outlierPtIds), f2) = zeros(numel(outlierPtIds), 1);

end

