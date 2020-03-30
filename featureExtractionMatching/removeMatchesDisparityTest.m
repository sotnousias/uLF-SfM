function [lfFeatures] = removeMatchesDisparityTest(lfFeatures)
% lfFeatures    -   NxM cell where N is the number of frames and M the
%                   number of features per frame
%
outlFeats = 0;

for n = 1:numel(lfFeatures) % for each frame
    for m = 1:numel(lfFeatures{n})
        
        xy = lfFeatures{n}{m}(:, 1:2);
        xy = xy';
        xy = xy(:);
        xy = xy';
        nproj=size(xy,2)/2; %25
        xy1=(reshape(xy, 2, nproj))'; % one projection per row
%         vproj=union(find(xy(:,1)~=-1), find(xy(:,2)~=-1)); % projections with -1 in x or y
        vxy=xy1;
        x=vxy(:,1);
        y=vxy(:,2);
        
        medx=median(x);
        madx=1.4826*median(abs(x-medx));
        robzscx=abs(x-medx)/madx;
        
        medy=median(y);
        mady=1.4826*median(abs(y-medy));
        robzscy=abs(y-medy)/mady;
        
        outl=find(max(robzscx, robzscy)>3.0);
        inl=find(max(robzscx, robzscy)<=3.0);
        
        lfFeatures{n}{m} = lfFeatures{n}{m}(inl, :);
        outlFeats = outlFeats + numel(outl);
%         xy(vproj(outl), :)=[-1 -1]; % remove outliers
%         xy=reshape(xy', 1, 2*nproj);

    end
end


fprintf('Found %d outliers\n', outlFeats);

end

