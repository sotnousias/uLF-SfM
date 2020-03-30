function mat = posesToMat(poses)
%

nposes = numel(poses);
mat = zeros(nposes, 7);

for n = 1:nposes
    R = poses{n}(1:3, 1:3);
    t = poses{n}(1:3, 4);
    t = t';
    q = rot2quat(R);
    mat(n, :) = [q' t];
    
end

end

