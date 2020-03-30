function [spatialRays] = lfPixelsToRays(pixelRays, intrinsicsLf, extrinsicsLf)
% map pixels in sub-aperture images to spatial rays using the intrinsics of
% the sub-apertures and the etrinsics of the sub-apertures.
% the input is a 4xM vector where the 1st two columns are the pixel and the
% last two are the x and y position of the sub-aperture in the grid
%
%
% intrinsicsLF -    3x3 intrisic matrix (a pinhole matrix)
%
% extrinsicsLF -    nxn extrinsic matrices obtained from LF calibtation
%                   (Yunsu Bok's toolbox for example)

%%%% TO DO: define the offset from the central sub-aperture image or input
%%%% the correct extrinsics matrix

numRays = size(pixelRays, 1);

spatialRays = zeros(size(pixelRays));

% get the pixels
pixels = pixelRays(:, 1:2)';

rays = pinv(intrinsicsLf) * [pixels; ones(1, size(pixels, 2))];
    
% 3rd coordinates already 1!
%rays(1:2, :) = rays(1:2, :) ./rays(3, :);
%%rays(1:2, :) = rays(1:2, :) ./repmat(rays(3, :), 2, 1);



for n = 1:numRays
   
    spatialRays(n, 1:2) = extrinsicsLf{pixelRays(n, 3), pixelRays(n, 4)}(1:2, 4)';
    spatialRays(n, 2) = -spatialRays(n ,2);
    spatialRays(n, 3:4) = rays(1:2, n)';
    
end




end

