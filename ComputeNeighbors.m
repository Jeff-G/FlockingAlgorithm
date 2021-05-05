% Computes a 2x? matrix, where each pair in the matrix represents neighboring indexes
% from the Nx3 points matrix.
% Inputs:
%   alg     : Enumerated NeighborAlgorithm value
%   points  : Nx3 matrix of 3D points
%   varargin: Additional algorithm-specific inputs
% Outputs:
%   neighbors: 2x? matrix of indexes into the points array that are neighbors
%   varargout: Additional algorithm-specific outputs
function [neighbors, varargout] = ComputeNeighbors(alg, points, varargin)
    if (alg == NeighborAlgorithm.Delaunay)
        % Compute the neighbors using Delaunay triangulation
        dt = delaunay(points);
        dtSorted = sort(dt, 2);
        neighbors = [dtSorted(:, 1) dtSorted(:, 2); ...
                     dtSorted(:, 1) dtSorted(:, 3); ...
                     dtSorted(:, 1) dtSorted(:, 4); ...
                     dtSorted(:, 2) dtSorted(:, 3); ...
                     dtSorted(:, 2) dtSorted(:, 4); ...
                     dtSorted(:, 3) dtSorted(:, 4)];
        neighbors = unique(neighbors, 'rows')';

        % Return the raw triangulation, if requested
        if (2 <= nargout)
            varargout{1} = dt;
        end
    elseif (alg == NeighborAlgorithm.Voronoi)
        % Compute the Voronoi diagram
        [vv, vc] = voronoin(points);

        % Find neighbors (there must be at least N edges)
        N = size(points, 1);
        neighbors = zeros(2, N);
        neighborCount = 0;
        for i = 1 : (N - 1)
            iIdxs = vc{i};
            for j = (i + 1) : N
                jIdxs = vc{j};
                % Compute the overlapping indexes
                oIdxs = intersect(iIdxs, jIdxs);
                % Uncomment this line to remove neighbors that only share one infinite face
                %oIdxs(oIdxs == 1) = [];
                % If the regions share a face, add the pair to the list of neighbors
                if (3 <= length(oIdxs))
                    neighborCount = neighborCount + 1;
                    neighbors(:, neighborCount) = [i; j];
                end
            end
        end

        % Return the raw tesselation results, if requested
        if (2 <= nargout)
            varargout{1} = vv;
            if (3 <= nargout)
                varargout{2} = vc;
            end
        end
    else
        error('Unsupported neighbor algorithm');
    end
end