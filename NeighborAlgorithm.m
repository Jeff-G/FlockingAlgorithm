classdef NeighborAlgorithm
    enumeration
        Delaunay % Delaunay triangulation neighbors
                 %   No additional inputs
                 %   Optionally outputs the result of calling the delaunay function
        Radius   % Forwards the call to Radius2 after squaring the optional input argument
        Radius2  % All points within the square of some radius
                 %   Optional argument for the square of the radius.
                 %     If not provided, all points are assumed to be neighbors
                 %   Optionally outputs a vector containing the square of the distance
                 %     between each pair of neighbors
        Voronoi  % Same as Delaunay, but much more computationally expensive.
                 %   No additional inputs
                 %   Optionally outputs the two results of calling the voronoin function
    end
end