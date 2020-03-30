function [raysOut] = rayCellToMatrix(raysIn)
% Convert rays which are in cell type to a matrix.


raysOut = [];

[dimx, dimy] = size(raysIn{1});



for n = 1:numel(raysIn)
    
    if (dimy == 6)
        raysOut = [raysOut; raysIn{n}];
    else
        
        raysOut = [raysOut raysIn{n}'];
    end
end
end

