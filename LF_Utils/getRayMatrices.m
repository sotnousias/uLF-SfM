function [raysL,raysR, rays1Ids, rays2Ids] = getRayMatrices(matches, sameSubview)
%UNTITLED18 Summary of this function goes here
%   Detailed explanation goes here

% sameSubview - boolean, set to 1 to use features from the same
%               sub-aperture images in the grid or set to 0 to use one to all, i.e. same
%               feature in frame i is in correspondence with all (matching) features of
%               frame j
%
%
% rays1Ids  -   Ids of each match, so outliers can be removed from the
%               matchesN struct.


raysL = [];
raysR = [];

idL = 0;
idR = 0;

if (sameSubview)
    idL = 1;
    idR = 1;
    
    for i = 1:numel(matches(1).rays3d)
        
        numRaysL = size(matches(1).rays3d{i}, 1);
        numRaysR = size(matches(2).rays3d{i}, 1);
        
        for j = 1:numRaysL
            
            % camera center in the grid, for the left frame 
            subCentL = matches(1).rays{i}(j, 3:4);
            
            % find the feature in the same sub-image
            dist = subCentL - matches(2).rays{i}(:, 3:4);
            
            dist = sqrt(sum(dist.^2, 2));
            
            idx = find(dist == 0);
            
            if (isempty(idx))
                continue;
            end
            
            raysL(idL, :) = matches(1).rays3d{i}(j, :);
            raysR(idR, :) = matches(2).rays3d{i}(idx, :);
            
            idL = idL + 1;
            idR = idR + 1;
        end
    end         
% use the one to all feature constrution    
else
    rays1Ids = [];
    rays2Ids = [];
    
    for i = 1:numel(matches(1).rays3d)
        
        numRaysL = size(matches(1).rays3d{i}, 1);
        numRaysR = size(matches(2).rays3d{i}, 1);
        
        for j = 1:numRaysL
            
            % each feature in the left light field, match it with all the
            % fearures from the right light field.
            raysL(idL+1:idL+numRaysR, :) = repmat(matches(1).rays3d{i}(j, :), numRaysR, 1);
            tmp = repmat([i j], numRaysR, 1);
            rays1Ids = [rays1Ids; tmp];
            
            raysR(idR+1:idR+numRaysR, :) = matches(2).rays3d{i};
            tmp = [i*ones(numRaysR, 1) (1:numRaysR)'];
            rays2Ids = [rays2Ids; tmp];
            idL = idL + numRaysR;
            idR = idR + numRaysR;
            
        end   
    end
end
end

