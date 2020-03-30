function [visMatrix, visBool, lfFeatures] = removeOutliersLfSfMAbsPose(visMatrix, visBool, ptsW, rays3, inlAbs, idsq, poses, frameId, lfFeatures, K )
% it looks like we have a problem with duplicates
% we also check reprojection error to remove outliers


ns = 1:size(rays3, 1);
outl = ns(~ismember(ns, inlAbs));

rm = idsq(outl, :);

Rh = poses{frameId}(1:3, 1:3);
th = poses{frameId}(1:3, 4);
ext = [Rh' -Rh'*th; zeros(1, 3) 1];

currId = 1;
outlReprTh = 1.5;

while (currId < size(rm, 1))
    
    % loop through all the points
    qids = find(rm(:, 1) == rm(currId, 1));
    
    pt3dId = rm(currId, 1);
%     if (pt3dId==1165)
%         5;
%     end
    
    % check if we already updated the visBool matrix and no visibility
    % exists for this correspondence
    if ~visBool(pt3dId, frameId)
        currId = currId + numel(qids) + 1;
        continue;
    end

    
    % pixels
    pixs =  lfFeatures{frameId}{visMatrix(pt3dId, frameId)}(rm(qids, 2), 1:2)';
    
    % reprojection errors
    reprs = zeros(numel(qids), 1);
    for n = 1:numel(qids)
        p = K * [eye(3) rays3(outl(qids(n)), 4:6)'] * ext * [ptsW(:, outl(qids(n))); 1];
        p = p ./ p(3);
        p = p(1:2);
        reprs(n) = norm(p - pixs(:, n));
    end
    
    % ids to be removed
    qs = find( reprs > outlReprTh);
    
    if (isempty(qs))
        currId = currId + numel(qids) + 1;

        continue;
    end
    
    % if the number of features to be removed are the total number of
    % features, just make zero the corresponding element in the visMatrix
    % and the same entry in visBool.
    if (numel(qs) == size(lfFeatures{frameId}{visMatrix(pt3dId, frameId)}, 1))
        visMatrix(pt3dId, frameId) = 0;
        visBool(pt3dId, frameId) = 0;
    else
        % make the outliers infinity so that we can keep the inliers in the
        % next step
        lfFeatures{frameId}{visMatrix(pt3dId, frameId)}(rm(qids(qs), 2), 1) = inf(numel(qs), 1);
        tmp = lfFeatures{frameId}{visMatrix(pt3dId, frameId)};
        tmp = tmp(~(tmp(:, 1) == inf), :);
        lfFeatures{frameId}{visMatrix(pt3dId, frameId)} = tmp;
        
        % now clear all the other presences of the same id since we have
        % this point as an inlier in the end (i.e. remove duplicates of
        % this featId)
        featId = visMatrix(pt3dId, frameId);
        
        prs = find(visMatrix(:, frameId) == featId);
        prs = prs(~ismember(prs, pt3dId));
        
        if isempty(prs)
            currId = currId + numel(qids) + 1;
            continue;
        else
            visMatrix(prs, frameId) = zeros(numel(prs), 1);
            visBool(prs, frameId) = zeros(numel(prs), 1);
            
        end
        
        % may be of use in debugging
        %     visualiseLfFeatures(lf{f1}, lfFeatures{f1}{visMatrix(w, f1)});
        %     visualiseLfFeatures(lf{nextFrame}, lfFeatures{nextFrame}{visMatrix(w, nextFrame)});
        %     visualiseLfFeatures(lf{nextFrame}, lfFeatures{nextFrame}{visMatrix(w, nextFrame)}(rm(qids(qs), 2), :), 'og')
        %     pause();
        
        currId = currId + numel(qids) + 1;
        
    end
    
end

