%#ok<*NASGU>
%#ok<*UNRCH>
function FinalProject()
    % Model Constants
    N       = 200; % Flock size
    deltaT  = 0.1; % Discrete model time step (s)
    finalT  = 120; % Simulation length (s)

    % Output settings
    traceN        = 3;     % Prints a trace for the first N points
    traceT        = inf;   % Maximum trace length (s)
    connectN      = 1;     % Connects the first N points to all of their neighbors
    identConvHull = true;  % Whether to uniquely draw birds on the convex hull
    genVideo      = true;  % Whether to generate a video of the resulting simulation
    genFinalFrame = true;  % Whether to save the final frame as an image
    genHullHist   = true;  % Whether to generate a histogram of the amount of time each
                           %     bird spends on the convex hull of the flock

    % Constants used to compute acceleration at each time step
    aMax     = 20; % Maximum acceleration (m/s^2)
    aIdeal   = 3;  % Acceleration toward vIdeal at vMin/vMax (m/s^2)
    aCenter  = 0.5;% Multiplied by the vector pointing to the center of all Voronoi neighbors
    dRepulse = 2;  % Distance beyond which the repulsive acceleration is zero (m)
    aRepulse = 1;  % Repulsion factor from nearby neighbors
    aRandom  = 1;  % Variance of random acceleration along each axis (scalar or 1x3 row vector)
    dAlign   = 4;  % Distance within which we begin trying to match alignment (m)
    aAlign   = 0.2;% Multiplied by the vector that would fully align a bird to neighbors

    % Altitude pressure
    hSpread = 20; % Distance from the ideal flight altitude (m) at which the acceleration
                  % toward hIdeal is 1 m/s^2

    % Velocity ranges for eastern crows are taken from the results in Table 1 from
    % "Flight-Speeds of Hawks and Crows", by Maurice Broun and Ben V. Goodwin, October 1943.
    vMin   =  8; % Minimum velocity for sustained flight (m/s, ~18 mph)
    vIdeal = 12; % Velocity that is easy to maintain over long distances (m/s, ~27 mph)
    vMax   = 16; % Maximum velocity (m/s, ~36 mph)
    vIdealVar = 0.25; % Ideal velocity variance (each bird has its own ideal velocity)

    % Lambda that computes the camera position at each time step
    azPeriod = 40; % Time to complete a 360 degree rotation (s)
    elPeriod = 30; % Time to cycle through an elevation range (s)
    getCameraViewAngles = @(i) deal( ...
        45 + (360 * deltaT / azPeriod) * i, ...                 % Azimuth
        15 + 20 * sin(deg2rad((360 * deltaT / elPeriod) * i))); % Elevation

    % Move the center of the flock slightly so velocity vectors are less likely to be cut off
    flockCenter = [-0.5 * deltaT * vIdeal, 0, 0];

    % Computed to generate a flock at some expected density
    pVariance = nthroot(N, 3); % Variance of the randomly generated points in p

    % Initialize the traces matrix (steps+1 x tracesN x 3)
    steps  = ceil(finalT / deltaT);
    traces = zeros(steps + 1, traceN, 3);

    % Compute the maximum number of trace points to display
    maxTraceLen = min(ceil(traceT / deltaT), steps);

    % N x steps+1 matrix that tracks whether a bird is on the convex hull of the flock
    isOnHull = false(N, steps + 1);

    % Compute constants used to calculate repulsive forces
    dRepulse2 = dRepulse^2;    % Square of the distance beyond which repulsion is zero
    oRepulse  = -1 / dRepulse; % Repulsive force added to the computed value, so the
                               % computed force is continuous over the dRepulse boundary

    % Compute constants used to calculate alignment forces
    dAlign2 = dAlign^2;    % Square of the distance beyond which no alignment is computed

    % Get the quadratic coefficients used to compute a force that tends to make all birds fly
    % at vIdeal.  This force is computed from a piecewise quadratic equation, so at vMin the
    % result is 1, at vIdeal it is zero, and at vMax it is -1.  This value is then multiplied
    % by aIdealMax to compute the resulting acceleration in the direction of the velocity.
    vMinCoef = [ 1 / (vIdeal - vMin)^2, 0, 0]; % If moving too slow, speed up
    vMaxCoef = [-1 / (vMax - vIdeal)^2, 0, 0]; % If moving too fast, slow down

    % Get the quadratic coefficients used to compute the force that pushes the flock toward
    % the ideal elevation
    hCoef = [1 / hSpread^2, 0, 0];

    % Compute the square of the maximum acceleration (aMax^2)
    aMax2 = aMax^2;

    % Generate a random 3d flock moving uniformly at vIdeal
    vIdeal = vIdeal + vIdealVar * randn(N, 1);
    p = pVariance .* randn(N, 3) + flockCenter;
    v = [vIdeal zeros(N, 2)];

    % Split the flock in half
    %No2 = floor(N / 2);
    %p(1:No2, 2) = p(1:No2, 2) + 40;

    % Get a string that includes the number of points and duration
    if (mod(finalT, 60) == 0)
        modelLen = sprintf('%dm', finalT / 60);
    else
        modelLen = sprintf('%ds', finalT);
    end
    modelStats = sprintf('%d@%s', N, modelLen);
    plotFileName = ['Birds Flocking (' modelStats ')'];

    % Create the axis and video
    if (genVideo || genFinalFrame)
        % Create the axis with fixed bounds
        fig = figure('Position', [0 0 1080 1080]);
        ax = gca;
        axis(ax, 8 * [-1, 1, -1, 1, -1, 1]);
        axis(ax, 'equal');
        ax.XLimMode = 'manual';
        ax.YLimMode = 'manual';
        ax.ZLimMode = 'manual';
        xlabel(ax, 'x');
        ylabel(ax, 'y');
        zlabel(ax, 'z');
        grid(ax, 'on');
        hold(ax, 'on');

        % Create a video of the flocking
        if (genVideo)
            video = VideoWriter([plotFileName '.avi']);
            video.FrameRate = 1 / deltaT;
            open(video);
        end
    end

    % Iterate over each time step
    for t = 0:steps
        % Compute the position at the next time step
        p = p + v * deltaT;

        % Center the flock at the desired flock center point
        p = p + (flockCenter - [mean(p(:, 1:2)) 0]);

        % Compute the current speed of each bird
        vMag = vecnorm(v, 2, 2);

        % Compute the normalized acceleration toward vIdeal
        vIdealDelta = vMag - vIdeal;
        tooSlow = vIdealDelta < 0;
        vIdealDelta( tooSlow) = polyval(vMinCoef, vIdealDelta( tooSlow));
        vIdealDelta(~tooSlow) = polyval(vMaxCoef, vIdealDelta(~tooSlow));
        vAccel = (aIdeal * vIdealDelta) .* (v ./ vMag);

        % Compute the altitude pressure
        hAccel = [zeros(N, 2) polyval(hCoef, p(:, 3)) .* -sign(p(:, 3))];

        % Compute the neighbors for each point
        neighbors = ComputeNeighbors(NeighborAlgorithm.Delaunay, p);

        % Compute an acceleration toward the center of visible flockmates
        nAccel = zeros(N, 3);
        for i = 1:N
            % Get a list of indexes into p that are neighbors to point i
            [row, col] = find(neighbors == i);
            attractNeighbors = neighbors(sub2ind(size(neighbors), 3 - row, col));

            % Compute the vector to the center of all neighbors
            centerVec = mean(p(attractNeighbors, :) - p(i, :), 1);
            nAccel(i, :) = aCenter * centerVec;
        end

        % Compute repulsion and alignment between flockmates
        aAccel = zeros(N, 3);
        rAccel = zeros(N, 3);
        for i = 1:N
            % Compute the square of the distance to each flockmate
            vecToMe = p(i, :) - p;
            distFromMe2 = sum(vecToMe .^ 2, 2);

            % Compute the flockmates that are within dAlign/dRepulse of this bird
            alignNeighbors = distFromMe2 < dAlign2;
            repulseNeighbors = distFromMe2 < dRepulse2;
            alignNeighbors(i) = false;
            repulseNeighbors(i) = false;

            % Compute the alignment forces
            if (any(alignNeighbors))
                % Compute the weighted average of the direction of nearby neighbors
                vNeighbors = v(alignNeighbors, :);
                vNeighbors = sum(vNeighbors .* (dAlign - ...
                    sqrt(distFromMe2(alignNeighbors))) ./ vecnorm(vNeighbors, 2, 2), 1);
                % Compute an acceleration that would result in fractional alignment with
                % these neighbors
                vNeighbors = vNeighbors * (norm(v(i, :)) / norm(vNeighbors));
                aAccel(i, :) = aAlign * (vNeighbors - v(i, :));
            end

            % Compute the repulsive forces
            if (any(repulseNeighbors))
                rForce  = vecToMe(repulseNeighbors, :) ./ distFromMe2(repulseNeighbors);
                rFactor = vecnorm(rForce, 2, 2);
                rForce  = rForce .* ((rFactor + oRepulse) ./ rFactor);
                rAccel(i, :) = aRepulse * sum(rForce, 1);
            end
        end

        % Compute the total acceleration (ensure it doesn't exceed the maximum allowed)
        a = vAccel + hAccel + nAccel + rAccel + aAccel + (aRandom .* randn(N, 3));
        aMag2 = sum(a .^ 2, 2);
        aTooHigh = aMax2 < aMag2;
        a(aTooHigh, :) = aMax * (a(aTooHigh, :) ./ sqrt(aMag2(aTooHigh)));

        % Compute the velocities at the next time step
        v = v + a * deltaT;

        % Compute the convex hull
        isOnHullArg = ComputeConvexHull(p);
        isOnHull(:, t + 1) = isOnHullArg;

        % Append the current point to the list of traces
        traceRange = (max(0, t - maxTraceLen):t) + 1;
        traces(traceRange(end), :, :) = p(1:traceN, :);

        % Generate a plot
        if (genVideo || (genFinalFrame && t == steps))
            % Clear the axis, and set the view angle
            cla(ax);
            [az, el] = getCameraViewAngles(t);
            view(ax, az, el);

            % Plot the points and velocities
            if (~identConvHull)
                isOnHullArg = [];
            end
            QuiverWithClassifier(ax, p, v * deltaT, isOnHullArg);

            % Draw the traces
            plot3(traces(traceRange, :, 1), ...
                  traces(traceRange, :, 2), ...
                  traces(traceRange, :, 3), '-');

            % Add neighbor connections
            for i = 1 : connectN
                ConnectNeighbors(ax, p, neighbors, i);
            end

            % Write the video frame
            drawnow;
            if (genVideo)
                writeVideo(video, getframe(fig));
            end
        end
    end

    % Restore normal plotting behavior
    if (genVideo || genFinalFrame)
        hold(ax, 'off');
    end

    % Save the final frame of the animation as an image
    if (genFinalFrame)
        TightenAxis(ax);
        saveas(fig, [plotFileName '.png'], 'png');
    end

    % Plot the fraction of time each bird spent on the convex hull
    if (genHullHist)
        % Create and label the plot
        histTitle = 'Flock Member Exposure to Predation';
        [fig, ax] = CreateAxis(histTitle);
        subtitle(ax, 'Amount of Time Each Bird Spent on the Flock''s Convex Hull', ...
            'FontSize', 16);
        xlabel(ax, sprintf('Time Exposed to Predation (%% of %s)', modelLen), ...
            'FontSize', 16, 'FontWeight', 'bold');
        ylabel(ax, sprintf('Birds (%d Total)', N), 'FontSize', 16, 'FontWeight', 'bold');

        % Plot the histogram of exposure to predation
        hullData = 100 * mean(isOnHull, 2);
        histHandle = histogram(ax, hullData);
        xMax = histHandle.BinWidth * ceil(max(hullData) / histHandle.BinWidth);
        xlim(ax, [0 xMax]);
        set(ax, 'XTick', 0:histHandle.BinWidth:xMax);

        % Compute and add the mean and standard deviation
        avgOnHull = mean(hullData);
        stdOnHull = std(hullData);
        yMax = ax.YLim(2);
        errorbar(ax, avgOnHull, 0.97 * yMax, stdOnHull, 'horizontal', 'o', ...
            'MarkerSize', 20, 'LineWidth', 2, 'CapSize', 20);
        extraText = { ...
            '\mu'   , sprintf(' = %.2f%%', avgOnHull); ...
            '\sigma', sprintf(' = %.2f%%', stdOnHull)};
        for i = 1:size(extraText, 1)
            text(ax, 0.01 * xMax, (1.03 - 0.04 * i) * yMax, extraText{i, 1}, ...
                'FontSize', 14, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
            text(ax, 0.02 * xMax, (1.03 - 0.04 * i) * yMax, extraText{i, 2}, ...
                'FontSize', 14, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
        end

        % Save the image
        TightenAxis(ax);
        saveas(fig, [histTitle ' (' modelStats ').png'], 'png');
    end
end