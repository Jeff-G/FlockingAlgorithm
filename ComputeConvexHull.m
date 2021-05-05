% Returns an Nx1 logical vector indicating if the associated point is on the convex hull
%   p: Nx3 matrix of points
function isOnHull = ComputeConvexHull(p)
    % Get the number of points
    N = size(p, 1);

    % Compute the convex hull
    hullTriangles = convhulln(p);

    % Identify the points on the convex hull
    hullPointIdxs = unique(reshape(hullTriangles, [], 1));
    isOnHull = false(N, 1);
    isOnHull(hullPointIdxs) = true;
end