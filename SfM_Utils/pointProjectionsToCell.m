function cl = pointProjectionsToCell(pts, pts3dIds, visMatrix, lfFeatures)

npts = numel(pts3dIds);
cl = cell(npts, 1);

LFsz = 1+25*2; % frame no + projections for one LF

for id = 1:npts
    ptId = pts3dIds(id);
    
    pt = pts(:, id)';
    
    lfSeen = find(visMatrix(ptId, :));
    
    numLfs = numel(lfSeen);
    projs=zeros(1, numLfs*LFsz);
    for i = 1:numLfs
        n = lfSeen(i);
        fs = -ones(25, 2);
        rs = lfFeatures{n}{visMatrix(ptId, n)}(:, 1:2);
%         fprintf('Point %d has %d projections in lf %d\n', id, size(rs, 1), n);
        subIds = lfFeatures{n}{visMatrix(ptId, n)}(:, 3:4);
        linIds = (subIds(:, 1) - 1 ) * 5 + subIds(:, 2);
        fs(linIds, :) = rs;
        fs = fs';
        fs = fs(:)';
        projs((i-1)*LFsz+1:i*LFsz)=[n-1 fs];
    end
    cl{id}=[pt numLfs projs];

end

end

