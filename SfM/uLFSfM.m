function [poses, pts, cols, pts3dIds, initIds, lfFeaturesOut, registered, visMatrix] = uLFSfM(lf, lfFeatures, visBool, visMatrix, geomVer, K, subLfExt, centDes, centFs, f1, f2, intrStr, lfIntrStr, reprTh, getDesc, frames, cameraSize, pyramidRes )
%

poses = cell(numel(lf), 1);
[visMatrix, visBool, lfFeatures, matchesN, pts3dIds, R, t] = initialiseMotion(lfFeatures, visBool, visMatrix, K, subLfExt, centDes, f1, f2);

initIds = 1;

poses{f1} = [eye(3) zeros(3, 1)];
poses{f2} = [R t];

pts3dIdsInit = pts3dIds;
matchesInit = matchesN;
lfFeaturesInitRel = lfFeatures;
visMatrixInitRel = visMatrix;



[pts3d, pts3D, ids, reprojectionErrs, pts3dIds, pts, cols, pc] = initialiseStructure(matchesN, poses, f1, f2, K, subLfExt, lf{f2}, reprTh, pts3dIds, 2, getDesc);


pts3dIdsInit = pts3dIds;
ptsInit = pts;
colsInit = cols;

registered = [f1 f2];
pts = ptsInit;
pts3dIds = pts3dIdsInit;
cols = colsInit;
ransacAbsTh = 0.25;
[width, height, ~] = size(lf{1}{3, 3});


queryImages = frames(~ismember(frames, registered));
notRegistered = [];
badFrames = [];
numRegistered = 0;
bundleBatch = 5;
looped = 0;

while (~isempty(queryImages))

    
    [pyramids, nextFrame] = selectNextView(pts3dIds, visMatrix, visBool, width, height, pyramidRes, queryImages, frames, centFs);
    
    if isnan(nextFrame) && ~looped
        queryImages = notRegistered;
        notRegistered = [];
        looped = 1;
    end
    
    if isnan(nextFrame) && looped
        fprintf('No more images to register...\n');
        break;
    end
    
    % find the ids of the already reconstructed 3d points that are seen in the
    % selected next frame
%     seenIds = find(visBool(pts3dIds, nextFrame));
    [Rw, tw, registerSucc, ptsW, rays3, idsq, inlAbs] = registerLfAbsolutePose(visBool, visMatrix, pts3dIds, nextFrame, pts, lfFeatures{nextFrame}, K, subLfExt);
    
    
    if ~(registerSucc)
        notRegistered = [notRegistered nextFrame];
        notRegistered = unique(notRegistered);
        fprintf('Could not register frame %d \n', nextFrame);
        queryImages = frames(~ismember(frames, registered));
        queryImages = queryImages(~ismember(queryImages, notRegistered));
        queryImages = queryImages(~ismember(queryImages, badFrames));

        continue;
    end
    poses{nextFrame} = [Rw tw];
    [visMatrix, visBool, lfFeatures] = removeOutliersLfSfMAbsPose(visMatrix, visBool, ptsW, rays3, inlAbs, idsq, poses, nextFrame, lfFeatures, K );
    
    
    % now try to reconstruct 3D points
    for regId = registered
        clear matchesN;
        
        qq = [regId nextFrame];
        qq = sort(qq);
        
        if ~(geomVer(qq(1), qq(2)))
%             fprintf('Skipping triangulation for unverified images..\n');
            continue;
        end
        [visBool, visMatrix, poses, pts, cols, pts3dIds, lfFeatures] = incrementalMapper(lf{nextFrame}, visBool, visMatrix, regId, nextFrame, lfFeatures, centDes, K, subLfExt, poses, pts3dIds, pts, cols, reprTh, registered );
       
    end
    
    [visMatrix, visBool, lfFeatures] = filterRegisteredFrameFeatures(visMatrix, visBool, pts, pts3dIds, poses, nextFrame, lfFeatures, K, subLfExt);

    registered = [registered nextFrame];
    
    if (numel(registered) == 3 ||( numel(registered) == 4))
        fprintf('Registered %d new images, calling global Bundle Adjustment...\n', numel(registered));
        [posesMat, ptsCellBA] = bundleAdjustmentHelper(pts, visMatrix, poses, lfFeatures, registered, notRegistered, pts3dIds);
        [poses, pts] = lfBundleAdjustmentMEX(intrStr, lfIntrStr, ptsCellBA, posesMat, poses, registered);
    end
    
    numRegistered = numRegistered + 1;
    notRegistered = notRegistered(~ismember(notRegistered, registered));
    queryImages = frames(~ismember(frames, registered));
    notRegistered = [];
    fprintf('Registered images: %d \n', numel(unique(registered)));
    if (numRegistered == bundleBatch)
        fprintf('Registered %d new images, calling global Bundle Adjustment...\n', bundleBatch);
        [posesMat, ptsCellBA] = bundleAdjustmentHelper(pts, visMatrix, poses, lfFeatures, registered, notRegistered, pts3dIds);
        [poses, pts] = lfBundleAdjustmentMEX(intrStr, lfIntrStr, ptsCellBA, posesMat, poses, registered);        
        numRegistered = 0;
    end
end

%%Final Global Bundle Adjustment

fprintf('Registered %d new images, calling global Bundle Adjustment...\n', numel(registered));
[posesMat, ptsCellBA] = bundleAdjustmentHelper(pts, visMatrix, poses, lfFeatures, registered, notRegistered, pts3dIds);
[poses, pts] = lfBundleAdjustmentMEX(intrStr, lfIntrStr, ptsCellBA, posesMat, poses, registered);
lfFeaturesOut = lfFeatures;



end

