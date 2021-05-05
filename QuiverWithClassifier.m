% Plots a scattergram of points with vectors
%   ax : Axis handle
%   p  : Nx3 matrix of locations
%   v  : Nx3 matrix of velocities
%   cat: Optional Nx1 vector that classifies the points into groups.
%        An empty vector or scalar is assumed to indicate no classification.
function QuiverWithClassifier(ax, p, v, cat)
    classes = unique(cat);
    if (len(classes) <= 1)
        % Plot all points without regard to classification
        quiver3(ax, p(:, 1), p(:, 2), p(:, 3), ...
            v(:, 1), v(:, 2), v(:, 3), 'AutoScale', false);
    else
        % Iterate over each class
        for class = classes
            % Get the points and velocities in this class
            isInClass = (cat == class);
            pInClass = p(isInClass, :);
            vInClass = v(isInClass, :);

            % Plot the points and velocities in this class
            quiver3(ax, pInClass(:, 1), pInClass(:, 2), pInClass(:, 3), ...
                vInClass(:, 1), vInClass(:, 2), vInClass(:, 3), 'AutoScale', false);
        end
    end
end