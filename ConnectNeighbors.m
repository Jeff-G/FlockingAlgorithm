% Draws lines to connect a point to all of its neighbors
%   ax       : Axis handle
%   points   : Nx3 matrix of points
%   neighbors: 2x? matrix of integer indices into the points array that are neighbors
%   idx      : Index into the points array whose neighbors will be drawn
function ConnectNeighbors(ax, points, neighbors, idx)
    for i = find(any(neighbors == idx))
        n = neighbors(:, i);
        line(ax, points(n, 1), points(n, 2), points(n, 3), 'Color', 'k');
    end
end