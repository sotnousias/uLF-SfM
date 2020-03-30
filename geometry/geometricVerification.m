function [geomVer] = geometricVerification(lfFeatures, lfDesc, essHomRatio, matchThresh, filterMatches, K)
%
if (nargin < 3)
    essHomRatio = 2;
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

centDes = cell(numel(lfFeatures), 1);
centFs = cell(numel(lfFeatures), 1);
geomVer = zeros(numel(lfFeatures));

for n = 1:numel(lfFeatures)
    
    tmp = zeros(128, numel(lfDesc{n}));
    tmp2 = zeros(2, numel(lfFeatures{n}));
    
    for m = 1:numel(lfDesc{n})
        tmp(:, m) = lfDesc{n}{m}(1, :)';
        tmp2(:, m) = lfFeatures{n}{m}(1, 1:2)';
    end
    
    centDes{n} = tmp;
    centFs{n} = tmp2;
end

frames = 1:numel(lfFeatures);

visitedFrames = [];
for f1 = 1:numel(lfFeatures)
    
    visitedFrames = [visitedFrames f1];
    
    for f2 = frames(~ismember(frames, visitedFrames))
        
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
            if (numel(matches1) < 100)
                match{f1}{f2} = [];
                continue;
            end
            
            [e, inle] = ransacFitEssentialHelper(matches1, matches2, K, K, ransacTh);
            [H, inlh] = ransacHomographyHelper(matches1, matches2, ransacTh*10^3);
            
            inlERat = numel(inle) / size(ids, 2);
            inlHRat = numel(inlh) / size(ids, 2);
            
            if (numel(inle)/numel(inlh) > essHomRatio)
                geomVer(f1, f2) = 1;
            end
            tmp = match{f1}{f2};
            match{f1}{f2} = zeros(2, numel(inle));
            match{f1}{f2}(1, :) = tmp(1, inle);
            match{f1}{f2}(2, :) = tmp(2, inle);
            
        end
    end
end
end

