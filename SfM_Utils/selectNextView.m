function [pyramid, nextFrame] = selectNextView(pts3dIds, visMatrix, visBool, width, height, pyramidRes, queryImages, frames, centFs)
% next best view selection based on score of the image as described in the
% COLMAP paper
%
% pts3dIds      -   Ids of the already trinagulated points in the visibilyt matrix indexing.
% visMatrix     -   Number of points by number of frames, consisting of the
%                   ids of the matches in each column.
% visBool       -   Boolean matric of the visibility per point per frame
%                   (usually visBool = visMatrix ~=0).
% width         -   Width of the image.
% height        -   Height of the image.
% pyramidRes    -   Resolution of the visibility pyramid.
% queryImages   -   The ids of the images from which the next view will be
%                   selected.
% frames        -   1:numberLfFrames
% centFs        -   The central keypoints, cell array where each cell holds
%                   the keypoints of each central sub-aperture image.

lvls = 2.^(1:pyramidRes) + 1;

% we can index each cell just by taking the divisiom between the step of
% the linspace and the corresponding keypoint index.
for i = 1:pyramidRes
   
    ws = linspace(1, width, lvls(i));
    hs = linspace(1, height, lvls(i));
    
    stepw(i) = ws(2) - ws(1);
    steph(i) = hs(2) - hs(1);
    
end

for n = 1:numel(frames) 
   for i = 1:pyramidRes
      
       pyramid{n}.cells{i} = zeros(lvls(i)-1);
       pyramid{n}.stepw(i) = stepw(i);
       pyramid{n}.steph(i) = steph(i);
       pyramid{n}.score(i) = lvls(i)-1;
       pyramid{n}.S = 0;
       pyramid{n}.frame = n;   
       
   end  
end


scores = zeros(numel(frames), 1);


for qid = queryImages
    
    seeIds = find(visBool(pts3dIds, qid));
    
    if isempty(seeIds)
        continue;
    end
    
    seeIds = pts3dIds(seeIds);
    
    
    fs = centFs{qid}(1:2, visMatrix(seeIds, qid));
    
    numFs = size(fs, 2);
    
    for n = 1:numFs
        for res = 1:pyramidRes
            c1 = floor(fs(1, n) / pyramid{qid}.steph(res)) + 1;
            c2 = floor(fs(2, n) / pyramid{qid}.stepw(res)) + 1;
            
            if (c1 > size(pyramid{qid}.cells{res}, 1))
                c1 = size(pyramid{qid}.cells{res}, 1);
            end
            
            if (c2 > size(pyramid{qid}.cells{res}, 1))
                c2 = size(pyramid{qid}.cells{res}, 1);
            end
            
            pyramid{qid}.cells{res}(c2, c1) = 1;
        end
    end
    
    
    for res = 1:pyramidRes
       pyramid{qid}.S = numel(find(pyramid{qid}.cells{res}))*pyramid{qid}.score(res) + pyramid{qid}.S;       
    end
    
    scores(qid) = pyramid{qid}.S;
    
end

[val, maxId] = max(scores);

if (val==0)
    nextFrame = NaN;
else

    nextFrame = pyramid{maxId}.frame;
end




end

