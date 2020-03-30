function [numMatches, visBool] = getPairwiseNumMatches(matches, visBool, badFrames)
%UNTITLED3 Summary of this function goes here

if (nargin < 3)
    badFrames = [];
end

numFrames = numel(matches) + 1;
numMatches = zeros(numFrames);

for i = 1:numFrames-1
    for j = 1:numFrames
        numMatches(i, j) = size(matches{i}{j}, 2);
%         numMatches(j, i) = numMatches(i, j);
    end
end

for i = 1:numFrames
    for j = 1:numFrames
        numMatches(j, i) = numMatches(i, j);
    end
end

if ~isempty(badFrames)
    for n = 1:numel(badFrames)
        numMatches(:, badFrames(n)) = zeros(numel(numMatches(:, badFrames(n))), 1);
        numMatches(badFrames(n), :) = zeros(1, numel(numMatches(badFrames(n), :)));
        visBool(:, badFrames(n)) = zeros(numel(visBool(:, badFrames(n))), 1);
    end
end


end

