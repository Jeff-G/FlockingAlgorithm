% Removes extraneous whitespace around a plot before saving it to disk
%   ax: Axis handle
function TightenAxis(ax)
    outerPos = ax.OuterPosition;
    tightIns = ax.TightInset;
    ax.Position = [...
        outerPos(1) + tightIns(1), ...
        outerPos(2) + tightIns(2), ...
        outerPos(3) - tightIns(1) - tightIns(3), ...
        outerPos(4) - tightIns(2) - tightIns(4), ...
    ];
end