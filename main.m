%% I. Get the pairwise matches
warning('off')

intrStr = './txtFiles/K.txt';
lfIntrStr = './txtFiles/subApertureRelPoses.txt';
% compute rootfist descriptors
getDesc = 1;

% number of levels in the pyramid of next view selection
pyramidRes = 3;
frames = 1:numel(lf);

% minimum number of views a feature is required to be seen in a LF frame
numValidViews = 4;

% threshold used for feature matching
matchThresh = 2.5;
sameSubviews = 0;
cameraSize = 1;
ransacTh = 1e-5;
essHomRatio = 3;
peakThresh = 0.0066;
onlyCentral = 0;
centId = 1;
badFrames = [];

% reprojection threshold
reprTh = 1.0;

lfFeatures = cell(numel(lf), 1);

for n = 1:numel(lf)

    lfFeatures{n} = extractLfFeaturesInFrame(lf{n}, numValidViews, matchThresh);
    
end

for n = 1:numel(lf)
[lfFeatures{n}, lfDesc{n}] = extractLfFeaturesInFrame(lf{n}, numValidViews);
end
%%

[~, matches, geomVer] =  sceneGraphV2(lfFeatures, lfDesc, essHomRatio, matchThresh, 1, K);

%% IV. Multi-image matches

[visMatrix, visBool] = getVisibilityMatrix(matches);

centFs = cell(numel(lf), 1);

[numMatches, visBool] = getPairwiseNumMatches(matches, visBool);

visBoolInit = visBool;
visMatrixInit = visMatrix;
lfFeaturesInit = lfFeatures;

[centFs, centDes] = getCentralFeatsDescrs(lfFeatures,  lfDesc);

[lfFeatures] = removeMatchesDisparityTest(lfFeatures);

fprintf('Finished feature extraction and matching...\n');

[f1, f2] = selectInitialFrames(geomVer, numMatches);
registered = [];
registered = [registered f1 f2];

%%
tic
[poses, pts, cols, pts3dIds, initIds, lfFeaturesOut, registered, visMatrixOut] = uLFSfM(lf, lfFeatures, visBool, visMatrix, geomVer, K, subLfExt, centDes, centFs, f1, f2, intrStr, lfIntrStr, reprTh, getDesc, frames, cameraSize, pyramidRes);
toc



