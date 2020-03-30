function [ptIntersect] = linePlaneIntersection(ptPlane, planeNormal, ptLine, lineDir)
%

% compute the distance along the line as in Wikipedia page

if ( size(ptPlane, 1) ~= 3)
    ptPlane = ptPlane';
end

if ( size(planeNormal, 1) ~= 3)
    planeNormal = planeNormal';
end

if ( size(ptLine, 1) ~= 3)
    ptLine = ptLine';
end

if ( size(lineDir, 1) ~= 3)
    lineDir = lineDir';
end


d = (ptPlane - ptLine)' * (planeNormal) / (lineDir' * planeNormal);


ptIntersect = ptLine + d * lineDir;


end

