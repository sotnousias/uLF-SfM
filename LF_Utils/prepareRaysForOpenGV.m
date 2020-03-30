function [raysOut] = prepareRaysForOpenGV(raysIn)
%UNTITLED19 Summary of this function goes here
% prepare rays for use in opengv
% 
% raysIn -  Nx4 rays in the uv-st parameterization. uv are the pinhole
%           center and the st are the direction 
%
% raysOut - Nx6 matrix. The first three elements of each row are the
%           direction of the ray, and the last three are the pinhole
%           center. Specifically:
%           - the st are padded with 1 in the end (still not sure if they
%           have to be normalised)
%           - the uv are padded with zero in the end
%
%

numRays = size(raysIn, 1);


% loop-free equivalent of the code below
%st = [raysIn(:, 3:4) ones(numRays,1)];
%st = st ./ repmat(sqrt(sum(st.^2,2)), 1, 3);
%raysOut = [st raysIn(:, 1:2) zeros(numRays,1)];


raysOut = zeros(numRays, 6);


for n = 1:numRays
    
    % pad the s-t coordinates
    raysOut(n, 1:2) = raysIn(n, 3:4);
    raysOut(n, 3) = 1;
    
    % normalize
    raysOut(n, 1:3) = raysOut(n, 1:3) ./ norm(raysOut(n, 1:3));
    
    % pad the u-v coordinates
    raysOut(n, 4:5) = raysIn(n, 1:2);
    raysOut(n, 6) = 0;
    
end


end

