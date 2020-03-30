function [visMatrix, matches, geomVer] = sceneGraphV2(lfFeatures, lfDesc, essHomRatio, matchThresh, filterMatches, K)
% OUTPUT
% V2 checks both f1 and f2 indices to see if new points will be appended
%
%
% matches     -     Cell array NxN where N is the the number of frames, each row
%                   corresponds to the matches between frame i and all the
%                   other frames. Each cell corresponds to the ids between
%                   frames i and j returned by the vl_ubcmatch. The ids
%                   correspond to the indices of the central
%                   features-descriptors.
%
% visMatrix   -     MxN matrix where M is the number of 3D points       
%                   (estimated in this function) and M is the number of frames.
%                   Each row contains the index of the central feature-descriptor.    

if (nargin < 3)
    essHomRatio = 3;
end

if (nargin < 4)
    matchThresh = 1.5;
end

if (nargin < 5)
    filterMatches = false;
    
else
    
    if (nargin < 6)
        error('K matrix is required for filtering matches with essential');
    else
        ransacTh = 1e-5;
        
    end
end

sigma = 0.8;

numframes=numel(lfFeatures);
centDes = cell(numframes, 1);
centFs = cell(numframes, 1);
geomVer = zeros(numframes);
match = cell(1, numframes-1);

for n = 1:numframes
    
    tmp = zeros(128, numel(lfDesc{n}));
    tmp2 = zeros(2, numel(lfFeatures{n}));
    
    for m = 1:numel(lfDesc{n})
        tmp(:, m) = lfDesc{n}{m}(1, :)';
        tmp2(:, m) = lfFeatures{n}{m}(1, 1:2)';
    end
    
    centDes{n} = tmp;
    centFs{n} = tmp2;
end

%%%%%% prepare for parfor
if mcrversion<=8
  if matlabpool('size') == 0 % check if pool is already open
     matlabpool('open', feature('numCores'));
  end
  %matlabpool('close');
else
  if isempty(gcp('nocreate')) % check if pool is already open
    parpool('local', 4);
  end
  %delete(gcp('nocreate'))
end

parfor f1 = 1:numframes-1
    
    for f2 = 1:numframes

        if (f2<=f1), continue; end % cannot have the above for start from f1+1 in outer parfor...
        
        % do symmetric matching?
        ids = vl_ubcmatch(centDes{f1}, centDes{f2}, matchThresh);
        ids2 = vl_ubcmatch(centDes{f2}, centDes{f1}, matchThresh);
        
        [~, ia, ib] = intersect(ids(2, :), ids2(1, :), 'stable');
        
        ids = ids(:, ia);
        match{f1}{f2} = ids;
        
        if (filterMatches)
            matches1 = centFs{f1}(1:2, match{f1}{f2}(1, :));
            matches2 = centFs{f2}(1:2, match{f1}{f2}(2, :));
            
            % check that you have enough points to get an estimate of the
            % essential matrix
            if (numel(matches1) < 50)
                match{f1}{f2} = [];
                continue;
            end
            
            [e, inle] = ransacFitEssentialHelper(matches1, matches2, K, K, ransacTh);
            [H, inlh] = ransacHomographyHelper(matches1, matches2, ransacTh*10^3);
            [F, inliF] = estimateFundamentalMatrix(matches1', matches2');
            [model, ~, ~] = gric_modsel(matches1, matches2, F, H, sigma);
            
            %if (numel(inle)/numel(inlh) > essHomRatio)
            %    geomVer(f1, f2) = 1;
            %end
            geomVer(f1, f2) = model;
            tmp = match{f1}{f2};
            match{f1}{f2} = zeros(2, numel(inle));
            match{f1}{f2}(1, :) = tmp(1, inle);
            match{f1}{f2}(2, :) = tmp(2, inle);
            
        end
    end

    %if (mod(f1, 10)==0), fprintf('sceneGraphV2: completed frame %d of %d\n', fr1, numframes); end
end

% corrs = [];
% corrs = zeros(size(match{1}{2}, 2), numel(lfFeatures));
% corrs (:, 1) = match{1}{2}(1, :)';
% corrs(:, 2) = match{1}{2}(2, :)';
% % ids1 = match{1}{2}(1, :);
% % ids2 = match{1}{2}(2, :);
% 
% for f1 = 1:numframes
%     for f2 = f1+1:numframes
%         
%         if (isempty(match{f1}{f2}))
%             continue;
%         end
%         
%         qids = match{f1}{f2};
%         newIds = 1:size(qids, 2);
%         [~, ia, ib] = intersect(qids(1, :), corrs(:, f1));
%         corrs(ib, f2) = qids(2, ia)';
%         newIds1 = newIds(~ismember(newIds, ia));
%         
%         [~, ia, ib] = intersect(qids(2, :), corrs(:, f2));
%         corrs(ib, f1) = qids(1, ia)';
%         newIds2 = newIds(~ismember(newIds, ia));
%         
%         newIds = intersect(newIds1, newIds2);
%         
%         tmp = zeros(numel(newIds), numel(lfFeatures));
%         tmp(:, f1) = qids(1, newIds);
%         tmp(:, f2) = qids(2, newIds);
%         corrs = [corrs; tmp];
%     end
% end

matches = match;
visMatrix = [];
geomVer(geomVer == 2) = 0;

end

