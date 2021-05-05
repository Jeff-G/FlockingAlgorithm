% Displays the Voronoi neighborhood around the point at the specified index
%   ax : Axis handle
%   vv : Voronoi vertices (?x3 list of 3D points)
%   vc : Cell array of indexes into the Voronoi vertices that comprise the neighborhood
%        for each point (Nx1 cell array of vectors containing integer indexes)
%   idx: Index into the vc cell array whose neighborhood will be displayed
function ShowVoronoi(ax, vv, vc, idx)
    % Remove the infinite vertex, if included in the list
    vvIdxs = vc{idx};
    if (vvIdxs(1) == 1)
        vvIdxs = vvIdxs(2:end);
    end

    % Display the convex hull of the remaining vertices
    faces = vv(vvIdxs, :);
    hull = convhulln(faces);
    trisurf(ax, hull, faces(:, 1), faces(:, 2), faces(:, 3), ...
        'FaceColor', [.8 .8 0], 'FaceAlpha', 0.2);
end