function [centFs, centDes] = getCentralFeatsDescrs(lfFeatures,  lfDesc)
%

centFs = cell(numel(lfFeatures), 1);

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


end

